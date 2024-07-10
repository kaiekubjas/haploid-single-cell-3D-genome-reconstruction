% This file shows some basic examples
% CVX is needed to use the function reconstruct_from_HiC_matrix 

% simulate a chromosome with 30 beads and generate the corresponding Hi-C matrix
radius = 2;
X0 = simulate_chromosomes(30);
A0 = generate_HiC_matrix(X0,r = radius);

% reconstruct the chromosome from the Hi-C matrix
[X1,flag] = reconstruct_from_HiC_matrix (A0, radius, dim = 3, method = 'unit_ball', local = "true",  backbone = "true", backbone_length=1);

% compute violation to original HiC-matrix A0
ob = objective_function(X1,A0,radius,method = 'unit_ball', backbone = 'false');
disp("objective function value is: "+string(ob));

% apply Procrustes transformation
[d,X] = procrustes(X0,X1,"scaling",false);
disp("Procrustes distance is: "+string(d));
