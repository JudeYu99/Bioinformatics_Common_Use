# Make the script executable by running the command:
#  chmod +x md_pca_cli.py
#
# Run the script with command-line arguments:
#  ./md_pca_cli.py md_data.npy --components pca_components.npy --variance explained_variance_ratio.npy


import argparse
from MDLovoFit.MDLovoFit import MDLovoFit
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='Perform PCA on MD data')
    parser.add_argument('md_data', type=str, help='Path to MD data file')
    parser.add_argument('--components', type=str, help='Path to save PCA components', default='pca_components.npy')
    parser.add_argument('--variance', type=str, help='Path to save explained variance ratio', default='explained_variance_ratio.npy')
    args = parser.parse_args()

    # load the MD data
    md_data = np.load(args.md_data)

    # create an instance of the MDLovoFit class
    md_pca = MDLovoFit(md_data)

    # perform PCA
    md_pca.perform_pca()

    # save the results
    np.save(args.components, md_pca.pca_components_)
    np.save(args.variance, md_pca.explained_variance_ratio_)

if __name__ == '__main__':
    main()
