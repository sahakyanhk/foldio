#!/usr/bin/env python3
"""Standalone ESMFold2 protein monomer structure prediction.

Two input modes (choose one):

  -ifasta FILE     A single FASTA file. Each record is predicted as a monomer
                   (no MSA). The record id (first token of the header) names
                   the output.

  -inputmsa FILE   A text file listing paths to MSA files (one path per line,
                   '#' comments allowed). Each MSA is an alignment in aligned
                   FASTA, a3m (.a3m), or Stockholm (.sto) format, chosen by
                   extension. Its first sequence is predicted as a monomer using
                   the whole alignment as the MSA. The MSA filename stem names
                   the output.

Output (into --output):
  <id>.cif      top predicted structure (mmCIF)
  metrics.tsv   id, length, plddt, ptm, n_msa  (also printed live)

Results stream one protein at a time and metrics are flushed as we go, so the
script handles large inputs and can be re-run to resume (finished .cif files
are skipped unless --overwrite).
"""

import sys
import argparse
from pathlib import Path

from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
from esm.models.esmfold2 import (
    ProteinInput,
    ESMFold2InputBuilder,
    StructurePredictionInput,
)
from esm.utils.msa.msa import MSA

GAPS = str.maketrans("", "", "-.")  # strip alignment gap characters


#==============================# INPUT PARSING #==============================#

def read_fasta(path: str):
    """Stream (id, sequence) records from a FASTA file, one at a time."""
    hid, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if hid is not None:
                    yield hid, "".join(seq)
                hid, seq = line[1:].split()[0], []
            else:
                seq.append(line)
    if hid is not None:
        yield hid, "".join(seq)


def fasta_queries(fasta_path: str):
    """Monomer queries from a plain FASTA: yields (id, sequence, None)."""
    for hid, seq in read_fasta(fasta_path):
        yield hid, seq.upper().translate(GAPS), None


def load_msa(msa_path: str) -> MSA:
    """Load an MSA, picking the parser from the file extension.

    .a3m -> a3m, .sto/.stockholm/.sth -> Stockholm, everything else (.fasta/
    .fa/.fas/.aln/.afa/...) -> aligned FASTA. The first record is the query.
    """
    ext = Path(msa_path).suffix.lower()
    if ext == ".a3m":
        return MSA.from_a3m(msa_path)
    if ext in (".sto", ".stockholm", ".sth", ".stk"):
        return MSA.from_stockholm(msa_path)
    aligned = [s.upper() for _, s in read_fasta(msa_path)]  # aligned FASTA, keep gaps
    return MSA.from_sequences(aligned)


def msa_queries(list_path: str):
    """One query per MSA file: yields (id, query_seq, MSA).

    Each listed path is an alignment (aligned FASTA / a3m / Stockholm); the
    first sequence is the query that gets folded.
    """
    for line in Path(list_path).read_text().splitlines():
        msa_path = line.strip()
        if not msa_path or msa_path.startswith("#"):
            continue
        msa = load_msa(msa_path)
        if msa.depth == 0:
            print(f"WARNING: empty MSA {msa_path}, skipping", file=sys.stderr)
            continue
        query = msa.query.upper().translate(GAPS)          # ungapped query to fold
        yield Path(msa_path).stem, query, msa


#=================================# PREDICT #=================================#

def rank_score(res) -> float:
    """Rank diffusion samples by pLDDT + pTM (both mapped to 0-1)."""
    plddt = float(res.plddt.mean())
    plddt = plddt / 100 if plddt > 1 else plddt
    ptm = float(res.ptm) if res.ptm is not None else 0.0
    return plddt + ptm


def predict_best(builder, model, seq, msa, args):
    """Fold one monomer and return the top-ranked diffusion sample."""
    prot = ProteinInput(id="A", sequence=seq, msa=msa)
    spi = StructurePredictionInput(sequences=[prot])
    out = builder.fold(
        model, spi,
        num_loops=args.num_loops,
        num_sampling_steps=args.num_steps,
        num_diffusion_samples=args.num_diffusion,
        seed=args.seed,
    )
    results = out if isinstance(out, list) else [out]
    return max(results, key=rank_score)


#===================================# MAIN #===================================#

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Standalone ESMFold2 protein monomer prediction",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("-ifasta", help="single FASTA file; each record folded as a monomer")
    g.add_argument("-inputmsa", help="file listing aligned-FASTA MSA paths (one per line)")
    p.add_argument("-output", required=True, help="output directory")
    p.add_argument("--num_diffusion", type=int, default=3,
                   help="diffusion samples per protein (best is kept)")
    p.add_argument("--num_loops", type=int, default=20, help="recycling loops")
    p.add_argument("--num_steps", type=int, default=200, help="sampling steps")
    p.add_argument("--seed", type=int, default=0, help="random seed")
    p.add_argument("--overwrite", action="store_true", help="re-fold existing outputs")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    queries = fasta_queries(args.ifasta) if args.ifasta else msa_queries(args.inputmsa)

    print("Loading ESMFold2 model ...", file=sys.stderr)
    model = ESMFold2Model.from_pretrained("biohub/ESMFold2").cuda().eval()
    builder = ESMFold2InputBuilder()

    summary = outdir / "metrics.tsv"
    write_header = not summary.exists()
    with open(summary, "a") as sf:
        if write_header:
            sf.write("id\tlength\tplddt\tptm\tn_msa\n")

        for n, (hid, seq, msa) in enumerate(queries, start=1):
            cif = outdir / f"{hid}.cif"
            if cif.exists() and not args.overwrite:
                print(f"[{n}] {hid}: exists, skipping")
                continue
            try:
                res = predict_best(builder, model, seq, msa, args)
            except Exception as e:
                print(f"[{n}] {hid}: FAILED ({e})", file=sys.stderr)
                continue

            cif.write_text(res.complex.to_mmcif())
            plddt = float(res.plddt.mean())
            ptm = float(res.ptm) if res.ptm is not None else float("nan")
            n_msa = msa.depth if msa is not None else 1

            sf.write(f"{hid}\t{len(seq)}\t{plddt:.4f}\t{ptm:.4f}\t{n_msa}\n")
            sf.flush()
            print(f"[{n}] {hid}  len={len(seq)}  plddt={plddt:.2f}  "
                  f"ptm={ptm:.3f}  n_msa={n_msa}  -> {cif.name}")

    print(f"Done. Metrics written to {summary}")


if __name__ == "__main__":
    main()
