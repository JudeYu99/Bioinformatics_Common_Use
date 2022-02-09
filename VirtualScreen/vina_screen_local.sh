#! /bin/bash

for f in fda*.pdbqt; do
    b=`basename $f .pdbqt`
    echo Processing ligand $b
    mkdir -p $b
    /home/toca/bioinfo/autodock/autodock_vina_1_1_2_linux_x86/bin/vina --config config.txt --ligand $f --out ${b}/out.pdbqt --log log_file/${b}_log.txt
done
