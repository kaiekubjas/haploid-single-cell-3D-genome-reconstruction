function [U] = generate_HiC_matrix(X,options)

% This function computes the HiC contact matrix according to our model
%
%
% Input:
% X            3-by-n matrix of coordinates of the chromosome X
% r (optional) radius in which an edge is observed 

arguments
    X
    options.r
end

if isfield(options,'r')
    r = options.r;
else
    r = 3*norm(X(1,:)-X(2,:));
end

n = size(X,1);

% Compute a matrix of contacts
U = zeros(n,n);
for i = 1:n
    for j = (i+1):n
        if norm(X(i,:)-X(j,:)) <= r
            U(i,j) = 1;
        else 
            U(i,j) = 0;
        end
        U(j,i) = U(i,j);
    end
end