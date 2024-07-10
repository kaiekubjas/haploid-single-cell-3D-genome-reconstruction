function [X] = simulate_chromosomes(n)

% This function simulates a chromose
%
% The chromosome consists of n beads, and the coordinates of these
% beads are recorded in a n-by-3 matrix called X.
%
% Input:
% n     desired number of beads


arguments
    n
end


% Simulate the chromosomes
X = zeros(n,3);
X(1,:) = [-1 -1 -1];
for i = 1:n-1
    step = randn(1,3);
    step = step/norm(step);
    X(i+1,:) = X(i,:)+step;
end

end
