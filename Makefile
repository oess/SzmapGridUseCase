SZMAP = szmap -mpi_np 8 -results_set max

default:
	@echo Please specify target.
	@cat Makefile | grep -E '^[A-Za-z0-9_.].*:' | grep -v '^default:' | awk -F: '{print $$1}' | sort | sed -e 's/\(.*\)/  -> \1/'

prep_all: split_schrodinger split_moe split_pch split_pch5
szmap_all: szmap_schrodinger szmap_moe szmap_pch szmap_pch5

all: prep_all szmap_all

split_schrodinger: data/5xco_schrodinger2017-4_minimized.mol2
	mkdir -p work
	python bin/split_complex.py -in data/5xco_schrodinger2017-4_minimized.mol2 \
		-p work/5xco_schrodinger2017-4_minimized_protein.oeb \
		-l work/5xco_schrodinger2017-4_minimized_ligand.oeb

split_moe: data/5xco_moe20160802_minimized.mol2
	mkdir -p work
	python bin/split_complex.py -in data/5xco_moe20160802_minimized.mol2 \
		-p work/5xco_moe20160802_minimized_protein.oeb \
		-l work/5xco_moe20160802_minimized_ligand.oeb

split_pch: data/5xco.pdb
	mkdir -p work
	mkhetdict data/5xco.pdb work/5xco_hets.txt
	reduce -db work/5xco_hets.txt -rotexist -build data/5xco.pdb > work/5xcoH.pdb 2> work/5xcoH_reduce.log || exit 0
	pch work/5xcoH.pdb work/5xcoH_complex.oeb work/5xcoH_others.oeb
	python bin/split_complex.py -in work/5xcoH_complex.oeb \
		-p work/5xco_pch_protein.oeb \
		-l work/5xco_pch_ligand.oeb

split_pch5: data/5xco.pdb
	mkdir -p work
	python bin/pch5.py -in data/5xco.pdb \
		-p work/5xco_pch5_protein.oeb \
		-l work/5xco_pch5_ligand.oeb

szmap_schrodinger: work/5xco_schrodinger2017-4_minimized_protein.oeb work/5xco_schrodinger2017-4_minimized_ligand.oeb
	$(SZMAP) \
		-prefix work/5xco_schrodinger \
		-protein work/5xco_schrodinger2017-4_minimized_protein.oeb \
		-around_mol work/5xco_schrodinger2017-4_minimized_ligand.oeb

szmap_moe: work/5xco_moe20160802_minimized_protein.oeb work/5xco_moe20160802_minimized_ligand.oeb
	$(SZMAP) \
		-prefix work/5xco_moe \
		-protein work/5xco_moe20160802_minimized_protein.oeb \
		-around_mol work/5xco_moe20160802_minimized_ligand.oeb

szmap_pch: work/5xco_pch_protein.oeb work/5xco_pch_ligand.oeb
	$(SZMAP) \
		-prefix work/5xco_pch \
		-protein work/5xco_pch_protein.oeb \
		-around_mol work/5xco_pch_ligand.oeb

szmap_pch5: work/5xco_pch5_protein.oeb work/5xco_pch5_ligand.oeb
	$(SZMAP) \
		-prefix work/5xco_pch5 \
		-protein work/5xco_pch5_protein.oeb \
		-around_mol work/5xco_pch5_ligand.oeb

###
### This is only for experiments

test_makegrid:
	if [ -z "${WEIGHT}" ] ; then \
		export WEIGHT=4 ; \
	fi ; \
	python experiments/merge.py -in work/5xco_pch5.oeb.gz -weight $${WEIGHT} -out shape.grd

test_rocs:
	python experiments/sim3D.py shape.grd testmc.oeb out.oeb

###


clean:
	rm -rf *~ *.bak .*~ .*.bak *.swp .*.swp *.pyc *.grd

complete_clean: clean
	rm -rf work
