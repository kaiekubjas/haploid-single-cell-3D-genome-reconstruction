# haploid-single-cell-3D-genome-reconstruction

This repository contains code for the article "Single cell 3D genome reconstruction in the haploid setting using rigidity theory" by Sean Dewar, Georg Grasegger, Kaie Kubjas, Fatemeh Mohammadi and Anthony Nixon.

## File descriptions

* MATLAB script `examples.m` that contains an example for synthetic data generation and 3D genome reconstruction from HiC matrix using our method.
* 
* MATLAB function `simulate_chromosomes.m` that outputs a synthetic chromosome of a given length in 3D. A chromosome of length n is represented by an (nx3)-matrix where each row represents coordinates of a bead on the chromosome.

* MATLAB function `generate_HiC_matrix.m` that outputs a HiC matrix given a chromosome and a threshold radius for a contact between beads.

* MATLAB function `reconstruct_from_HiC_matrix.m` that reconstructs a chromosome from a HiC matrix.

* MATLAB function `objective_function.m` that constructs the objective function for improving the reconstruction using gradient descent. The same function is used to construct the measure of dissimilarity to evalute the goodness of the reconstruction.

* Mathematica script `RealDataGraphGeneration.wl` to preprocess real data to obtain the HiC matrix in the format that can be used with our function `reconstruct_from_HiC_matrix.m`.

## Dependencies

The MATLAB script `examples.m` and function `reconstruct_from_HiC_matrix.m` require CVX: Matlab Software for Disciplined Convex Programming.
