# modelcraft-paper

## data/af

Contains 3030 AlphaFold MR structures.
The target structures only have a single protein entity,
although there could be multiple copies.
Each folder contains:

- contents.json
- data.mtz
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

The MTZ contains:

- FP,SIGFP
- FREE
- HLA,HLB,HLC,HLD (deposited experimental phases without density modification)
