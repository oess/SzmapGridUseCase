# SzmapGridUseCase
Explore the use case of Szmap grid

### Data

| filename | Hydrogens | Partial charges | Water |
|:---|:---:|:---:|:---:|
| data/5xco_schrodinger2017-4.mol2           | Y | N | N |
| data/5xco_schrodinger2017-4_minimized.mol2 | Y | Y | N |
| data/5xco_moe20160802.mol2                 | Y | Y | Y |
| data/5xco_moe20160802_minimized.mol2       | Y | Y | Y |
| data/5xco.pdb                              | N | N | Y |

### Preparation for Szmap
1) Water is available (forcibly eliminated)
2) Partial charges are given (if not available)
3) Hydrogens are placed (if not available)

See Makefile for prep in each individual case.

### Szmap execution
```
szmap -mpi_np 8 -results_set max -prefix work/PREFIX -protein PROTEINFILE.oeb -around_mol LIGANDFILE.oeb
```
