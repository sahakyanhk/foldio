# foldio

Misc scripts for structure prediction and search


install esmfold2, requares python=3.12
```bash
conda create -n esmfold2 python=2.12
pip install esm@git+https://github.com/Biohub/esm.git@main

#install search tools
bash install.sh 
```

## run esmfold 
Locally
``` 
#from a fasta file with multiple sequences
python structure_prediction/esmfold2.py -ifasta examples/proteins.fasta -o examples/proteins_output

#from a list of paths pointing to MSAs, structure for the fist sequence from each MSA is predicted
python structure_prediction/esmfold2.py -inputmsa examples/msas.list -o examples/msas_output
```

Submit an sbatch job



## search structures with foldseek

```bash
structure_search/foldseek_search.sbatch examples/proteins_output examples/foldseek_search_afdb.tsv
```

