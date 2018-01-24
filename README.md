# SzmapGridUseCase
Explore the use case of Szmap grid

### Data
(H: Hydrogens, C: Partial charges, W: Water)<br/>
data/5xco_schrodinger2017-4.mol2           H+, C-, W-<br/>
data/5xco_schrodinger2017-4_minimized.mol2 H+, C+, W-<br/>
data/5xco_moe20160802.mol2                 H+, C+, W+<br/>
data/5xco_moe20160802_minimized.mol2       H+, C+, W+<br/>
data/5xco.pdb                              H-, C-, W+<br/>

### Preparation for Szmap
1) Water is available (eliminated if)
2) Partial charges are given (if not available)
3) Hydrogens are placed (if not available)

### Szmap execution
```
szmap -mpi_np 8 -results_set max -prefix work/PREFIX -protein PROTEINFILE.oeb -around_mol LIGANDFILE.oeb
```
