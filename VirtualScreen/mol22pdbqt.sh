#! /usr/bin

for i in {1..2115}
do      /home/toca/bioinfo/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py  -l fda$i.mol2 -o fda$i.pdbqt
	echo $i
done
