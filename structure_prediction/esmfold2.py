#!/usr/bin/env python3
"""Standalone ESMFold2 structure prediction.

Input modes (choose exactly one):

  -ifasta FILE     A single FASTA file. Each record is folded as a protein
                   monomer (no MSA). The record id (first header token) names
                   the output.

  -inputmsa FILE   A text file listing paths to MSA files (one per line, '#'
                   comments allowed). Each MSA is an alignment in aligned FASTA,
                   a3m (.a3m), or Stockholm (.sto) format, chosen by extension.
                   Its first sequence is folded as a monomer using the whole
                   alignment as the MSA. The MSA filename stem names the output.

  -json FILE       A JSON file describing one or more queries, each of which may
                   be a single chain or a complex of protein / dna / rna / ligand
                   chains. Format:

                     {"queries": {
                        "<name>": {"chains": [
                          {"molecule_type": "protein",
                           "chain_ids": ["A", "B"],        # copies -> homomer
                           "sequence": "MKL..."},
                          {"molecule_type": "ligand",
                           "chain_ids": ["F", "G"],
                           "ccd_codes": "ATP"},            # or ["ATP", ...]
                          {"molecule_type": "ligand",
                           "chain_ids": "Z",
                           "smiles": "CC(=O)..."}
                        ]}}}

                   chain_ids may be a single id or a list; a list of N ids means
                   N identical copies. Each query <name> names the output.

Output (into -output):
  <name>.cif    top predicted structure (mmCIF)
  metrics.tsv   name, n_chains, plddt, ptm, iptm  (also printed live)

Queries stream one at a time and metrics are flushed as we go, so the script
handles large inputs and can be re-run to resume (finished .cif files are
skipped unless --overwrite).
"""

import sys
import json
import argparse
from pathlib import Path

from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
from esm.models.esmfold2 import (
    ProteinInput,
    RNAInput,
    DNAInput,
    LigandInput,
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
    """Yield (name, [ProteinInput]) monomers from a plain FASTA."""
    for hid, seq in read_fasta(fasta_path):
        prot = ProteinInput(id="A", sequence=seq.upper().translate(GAPS))
        yield hid, [prot]


def load_msa(msa_path: str) -> MSA:
    """Load an MSA, picking the parser from the file extension.

    .a3m -> a3m, .sto/.stockholm/.sth/.stk -> Stockholm, everything else
    (.fasta/.fa/.aln/.afa/...) -> aligned FASTA. The first record is the query.
    """
    ext = Path(msa_path).suffix.lower()
    if ext == ".a3m":
        return MSA.from_a3m(msa_path)
    if ext in (".sto", ".stockholm", ".sth", ".stk"):
        return MSA.from_stockholm(msa_path)
    aligned = [s.upper() for _, s in read_fasta(msa_path)]  # aligned FASTA, keep gaps
    return MSA.from_sequences(aligned)


def msa_queries(list_path: str):
    """Yield (name, [ProteinInput]) monomers, one per listed MSA file."""
    for line in Path(list_path).read_text().splitlines():
        msa_path = line.strip()
        if not msa_path or msa_path.startswith("#"):
            continue
        msa = load_msa(msa_path)
        if msa.depth == 0:
            print(f"WARNING: empty MSA {msa_path}, skipping", file=sys.stderr)
            continue
        query = msa.query.upper().translate(GAPS)          # ungapped query to fold
        prot = ProteinInput(id="A", sequence=query, msa=msa)
        yield Path(msa_path).stem, [prot]


def build_chain(chain: dict):
    """Turn one JSON chain spec into an ESMFold2 input object."""
    mtype = chain["molecule_type"].lower()
    ids = chain["chain_ids"]                       # str or list -> copies

    if mtype == "protein":
        return ProteinInput(id=ids, sequence=chain["sequence"].upper())
    if mtype == "dna":
        return DNAInput(id=ids, sequence=chain["sequence"].upper())
    if mtype == "rna":
        return RNAInput(id=ids, sequence=chain["sequence"].upper())
    if mtype == "ligand":
        if chain.get("smiles"):
            return LigandInput(id=ids, smiles=chain["smiles"])
        ccd = chain.get("ccd_codes") or chain.get("ccd")
        if not ccd:
            raise ValueError("ligand chain needs 'smiles' or 'ccd_codes'")
        return LigandInput(id=ids, ccd=[ccd] if isinstance(ccd, str) else list(ccd))
    raise ValueError(f"unknown molecule_type: {mtype!r}")


def json_queries(json_path: str):
    """Yield (name, [input objects]) for each query in a JSON file."""
    data = json.loads(Path(json_path).read_text())
    for name, spec in data["queries"].items():
        yield name, [build_chain(c) for c in spec["chains"]]


def count_chains(chains) -> int:
    """Total chain count, expanding chain_ids lists (copies)."""
    return sum(len(c.id) if isinstance(c.id, list) else 1 for c in chains)


#=================================# PREDICT #=================================#

def rank_score(res) -> float:
    """Rank diffusion samples by pLDDT + pTM + ipTM (all mapped to 0-1)."""
    plddt = float(res.plddt.mean())
    plddt = plddt / 100 if plddt > 1 else plddt
    ptm = float(res.ptm) if res.ptm is not None else 0.0
    iptm = float(res.iptm) if res.iptm is not None else 0.0
    return plddt + ptm + iptm


def predict_best(builder, model, chains, args):
    """Fold one query (monomer or complex); return top-ranked diffusion sample."""
    spi = StructurePredictionInput(sequences=chains)
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
        description="Standalone ESMFold2 structure prediction",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("-ifasta", "-i", dest="ifasta",
                   help="single FASTA file; each record folded as a monomer")
    g.add_argument("-inputmsa", "-m", dest="inputmsa",
                   help="file listing MSA paths (one per line); fold 1st seq of each")
    g.add_argument("-json", "-j", dest="json",
                   help="JSON file of queries (monomers/complexes, multi-molecule)")
    p.add_argument("-output", "-o", dest="output", required=True,
                   help="output directory")
    p.add_argument("--num_diffusion", type=int, default=3,
                   help="diffusion samples per query (best is kept)")
    p.add_argument("--num_loops", type=int, default=10, help="recycling loops")
    p.add_argument("--num_steps", type=int, default=100, help="sampling steps")
    p.add_argument("--seed", type=int, default=0, help="random seed")
    p.add_argument("--overwrite", action="store_true", help="re-fold existing outputs")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.ifasta:
        queries = fasta_queries(args.ifasta)
    elif args.inputmsa:
        queries = msa_queries(args.inputmsa)
    else:
        queries = json_queries(args.json)

    print("Loading ESMFold2 model ...", file=sys.stderr)
    model = ESMFold2Model.from_pretrained("biohub/ESMFold2").cuda().eval()
    builder = ESMFold2InputBuilder()

    summary = outdir / "metrics.tsv"
    write_header = not summary.exists()
    with open(summary, "a") as sf:
        if write_header:
            sf.write("name\tn_chains\tplddt\tptm\tiptm\n")

        for n, (name, chains) in enumerate(queries, start=1):
            cif = outdir / f"{name}.cif"
            if cif.exists() and not args.overwrite:
                print(f"[{n}] {name}: exists, skipping")
                continue
            try:
                res = predict_best(builder, model, chains, args)
            except Exception as e:
                print(f"[{n}] {name}: FAILED ({e})", file=sys.stderr)
                continue

            cif.write_text(res.complex.to_mmcif())
            nc = count_chains(chains)
            plddt = float(res.plddt.mean())
            ptm = float(res.ptm) if res.ptm is not None else float("nan")
            iptm = float(res.iptm) if res.iptm is not None else float("nan")

            sf.write(f"{name}\t{nc}\t{plddt:.4f}\t{ptm:.4f}\t{iptm:.4f}\n")
            sf.flush()
            print(f"[{n}] {name}  chains={nc}  plddt={plddt:.2f}  "
                  f"ptm={ptm:.3f}  iptm={iptm:.3f}  -> {cif.name}")

    print(f"Done. Metrics written to {summary}")


if __name__ == "__main__":
    main()
