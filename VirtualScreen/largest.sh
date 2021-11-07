export script_path=/Users/YuZhu99/Downloads/mgltools_i86Darwin9_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24
export python_path=/Users/YuZhu99/Desktop/

for name in `cat mol.lst`
do
#remove .mol2 from the name
  name=`basename ${name} .mol2`
  cd ${name}
 ${python_path}/pythonsh ${script_path}/write_largest_cluster_ligand.py -o ${name}_BC.pdbqt -v
  cd ..
done