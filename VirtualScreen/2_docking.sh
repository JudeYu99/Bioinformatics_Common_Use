#! /usr/bin/env bash

#_____________________readme first carefully !________________________________
#file name   : 02_docking.sh
#file needed : 
#              (1) receptor file in pdbqt format and receptor maps.
#              (2) many separate small molecules in mol2 format and mol.lst contain the name of these small molecules.
#
#Usage       :  02_docking.sh  -r [receptor.pdbqt] -l [mol list] -c [ cpu number that you use ] [-g num num_of_GA_run ] -p your_uniq_name
#
#Author      : zjxu@mail.shcnc.ac.cn
#Created     : 2/23/2012
#Last updated: 2/23/2012
#Last updated: 3/12/2012,                echo "Tips: If you are using ADT, I recommend you set map types directly: 'A C HD N NA OA S SA F Cl CL Br BR I'. Receptor grid maps were needed, either prepared from 01_grid.sh or from ADT by yourslef."


#_____________________________________________________________________________

USAGE()
{
        echo "-------------------------------------------------------"
        echo "./02_docking.sh -r [receptor.pdbqt] -l [mol list] -c [ cpu number that you use ] [-g num num_of_GA_run ] -p your_uniq_name"
        echo ""
        echo "for example:"
        echo "./02_docking.sh -r receptor.pdbqt -l mol.lst -c 4 [-g 20] -p zjxu"
        echo ""
        echo "Tips: If you are using ADT, I recommend you set map types directly: 'A C HD N NA OA S SA F Cl CL Br BR I'
Receptor grid maps were needed, either prepared from 01_grid.sh or from ADT by yourslef."

        echo "-------------------------------------------------------"

}


#change this to your own path
export script_path=/Library/MGLTools/1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24
export python_path=/Library/MGLTools/1.5.6/bin

num_of_GA_run=10
program_state=0

while getopts r:l:c:g:p: OPT
do
        case $OPT in
                r) receptor=${OPTARG%.pdbqt}
                   echo "receptor is ${receptor}.pdbqt"
                   program_state=1
                   ;;
                l) mol_file=${OPTARG}
                   echo "ligand list file is ${mol_file}"
                   ;;
                c) num_of_job=${OPTARG}
                   echo "number of parallel jobs is ${num_of_job}"
                   ;;
                g) num_of_GA_run=${OPTARG}
                   echo "number of GA runs is ${num_of_GA_run}"
                   ;;
                p) YOUR_UNIQ_NAME=${OPTARG}
                   program=`which autodock4`
                   cp $program ./autodock4_${YOUR_UNIQ_NAME}
                   ;;
                \?) USAGE
                   exit
                   ;;
        esac
done
if [ $program_state -eq 0 ]
then USAGE; exit 1
fi

#--------------------------------------------------------------------------
# do some tests

if [ ! -f ${receptor}.pdbqt ]
 then
  echo "can't find file ${receptor}.pdbqt"
  USAGE
  exit
fi

if [ ! -f ${mol_file} ]
 then
  echo "can't find file ${mol_file}"
  USAGE
  exit
fi

#--------------------------------------------------------------------------



#ligad preparation
for name in `cat mol.lst`
do
#remove .mol2 from the name
# name=`basename ${name} .mol2`
  name=${name%.mol2}
#generate ligand pdbqt file
#-d means dictionary
  ${python_path}/pythonsh ${script_path}/prepare_ligand4.py -l ${name}.mol2 -d ligand_dict.py

# Create a subdirectory named <ligand> and populate
# it with the docking input files: a) the pdbqt from the
# Ligands directory will be copied directly; b) the maps
# will be linked to the Receptor directory; and, c) the dpf
# file will be created using prepare_dpf4.py:
  if [ -e ${name}.pdbqt ]
    then
      mkdir ${name}
      cd ${name}
      cp ../${name}.pdbqt .
      ln -s ../${receptor}.pdbqt .
      ln -s ../${receptor}*map* .
      # generate dpf for each ligand
      ${python_path}/pythonsh ${script_path}/prepare_dpf42.py -l ${name}.pdbqt -r ${receptor}.pdbqt -p ga_num_evals=1750000 -p ga_pop_size=150 -p ga_run=${num_of_GA_run} -p rmstol=2.0

      # get the number of current running autodock4
      current_job_num=`ps -f -U ${USER} | grep -v 'grep' | grep -c "autodock4_${YOUR_UNIQ_NAME}"`

      # wait until ${current_job_num} is less than ${num_of_job}
      while [ ${current_job_num} -ge ${num_of_job} ]
        do
          sleep 10
          current_job_num=`ps -f -U ${USER} | grep -v 'grep' | grep -c "autodock4_${YOUR_UNIQ_NAME}"`
      done

      # run docking for each ligand
      #autodock4 -p ${name}_receptor.dpf -l ${name}.dlg
      ../autodock4_${YOUR_UNIQ_NAME} -p ${name}_${receptor}.dpf -l ${name}.dlg &
      cd ..
  fi
done

rm autodock4_${YOUR_UNIQ_NAME}
#To view the liand info, you should run ${script_path}/examine_ligand_dict.py > summary.txt
