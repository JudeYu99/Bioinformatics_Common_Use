#! /usr/bin/env bash
#change this to your own path
export script_path=/Users/YuZhu99/Downloads/mgltools_i86Darwin9_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24
export python_path=/Users/YuZhu99/Desktop/


for name in `cat mol.lst`
do
#remove .mol2 from the name
  name=`basename ${name} .mol2`
# Extract the Free Energy of Binding for the lowest energy
# in the largest cluster from the dlg files using the python
# script summarize_results4.py :
  ${python_path}/pythonsh ${script_path}/summarize_results4.py -d ${name} -t 2.0 -L -a -o summary_2.0.txt
  ${python_path}/pythonsh ${script_path}/write_clustering_histogram_postscript.py -d ${name} -t 8 -o ${name}_cluster.ps -v

#(1) read all the files with extension '.dlg' into one Docking
#(2) compute a clustering at the specified rms tolerance 
#(3) write the ligand with the coordinates of the lowest-energy conformation in the largest cluster to a file
done

#convert ps to png
for name in `ls *.ps`
do
  name=${name%.ps}
  convert ${name}.ps ${name}.png
done

#-k5n means sort on field 5 here the lowest energy in the largest cluster.  -t, means to use commas as field separators
cat summary_2.0.txt | sort -k5n -t, > summary_2.0.sort