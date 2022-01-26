# modelcraft-paper

## Software versions

CCP4 version 7.1.018

Python 3.8.10 with PIP packages installed:

- gemmi==0.5.1
- matplotlib==3.3.3
- modelcraft==2.2.3
- numpy==1.19.5
- pandas==1.0.4
- requests==2.22.0
- solrq==1.1.1

## data/af

Contains 3030 AlphaFold MR structures.
The target structures only have a single protein entity,
although there could be multiple copies.
Each folder contains:

- contents.json
- data.mtz
- deposited.cif.gz
- model.cif

The model is the AlphaFold model
where all the copies of the chains have been placed using MOLREP.
PLDDT cutoffs of 0, 50, 70 and 90 were all tried and the resulting structure
is the one with the highest MR score (so long as all copies could be placed).
This was done on the whole AlphaFold DB for the human proteome,
but only for models that weren't split between multiple files.

The MTZ contains:

- Observations (Ianom, Imean, Fanom or Fmean)
- FreeR_flag

## data/ep

Contains 1204 JCSG structures.
Each folder contains:

- contents.json
- data.mtz
- deposited.cif.gz

The MTZ contains:

- FP,SIGFP
- FREE
- HLA,HLB,HLC,HLD (deposited experimental phases without density modification)
