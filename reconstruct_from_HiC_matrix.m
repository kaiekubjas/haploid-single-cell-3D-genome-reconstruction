function [X,flag] = reconstruct_from_HiC_matrix(A,r,options)

% Function that computes the 3D reconstruction from a Hi-C matrix using
% semidefinite programming 
% 
% CVX is needed to use this function
%
% Input:
% A      Hi-C matrix
% r      radius in which an edge is observed 
% dim    reconstruction dimension, default 3
% method reconstruction method, options 'unit_ball', 'equalities_and_inequalities', 'marble_graph'

arguments
    A
    r
    options.dim
    options.method
    options.backbone
    options.backbone_length
    options.edges
    options.edge_lengths
    options.lambda
    options.epsilon
    options.local
end

if isfield(options,'dim')
    dim = options.dim;
else
    dim = 3;
end

if isfield(options,'method')
    method = options.method;
else
    method = 'unit_ball';
end

if isfield(options,'edges')
    edges = options.edges;
else
    edges = [];
end

if isfield(options,'edge_lengths')
    edge_lengths = options.edge_lengths;
else
    edge_lengths = [];
end

if isfield(options,'lambda')
    lambda = options.lambda;
else
    lambda = 0;
end

if isfield(options,'backbone')
    backbone = options.backbone;
else
    backbone = 'false';
end

if isfield(options,'backbone_length')
    backbone_length = options.backbone_length;
else
    backbone_length = 1;
end

% reset edges for methods which should not use them
if strcmp(method,'marble_graph') || strcmp(method,'unit_ball')
    edges = [];
end

if isfield(options,'epsilon')
    epsilon = options.epsilon;
else
    epsilon = 0.001;
end

if isfield(options,'local')
    local = options.local;
else
    local = 'false';
end

n = size(A,1);

% construct backbone edges
if backbone == "true"
    back_edges = zeros(n-1,2);
    for i = 2:n
        edge_lengths(i-1,i) = backbone_length^2;
        edge_lengths(i,i-1) = backbone_length^2;
        back_edges(i-1,1:2) = [i-1,i];
    end
else
    back_edges=[];
end

marble_edges = [];    
if strcmp(method,'marble_graph')
    for i = 1:n
        for j = 1:(i-1)
            if A(i,j) == 1
                edge_lengths(i,j) = r^2;
                edge_lengths(j,i) = r^2; 
                marble_edges = [marble_edges; [i,j]];                
            end
        end
    end
end


% make sure that edges are not treated as inequalities
A_edges = A;
for k = 1:size(edges,1)
    i = edges(k,1);
    j = edges(k,2);
    A_edges(i,j) = -1;
    A_edges(j,i) = -1;
end
A_back = A_edges;
for k = 1:size(back_edges,1)
    i = back_edges(k,1);
    j = back_edges(k,2);
    A_back(i,j) = -2;
    A_back(j,i) = -2;
end
A_marble = A_back;
for k = 1:size(marble_edges,1)
    i = marble_edges(k,1);
    j = marble_edges(k,2);
    A_marble(i,j) = -3;
    A_marble(j,i) = -3;
end
A_local = A_marble;
edges_input = edges;
edges_pure = [edges;back_edges];
edges = [edges_pure;marble_edges];

% main part
cvx_begin
    variable G(n,n) symmetric
    obj = -lambda*trace(G);
    for k = 1:size(edges,1)
        i = edges(k,1);
        j = edges(k,2);
        obj = obj + (G(i,i)+G(j,j)-2*G(i,j) - edge_lengths(i,j))^2; % note: edge_lengths are actually squared edge lengths
    end
    minimize obj;
    subject to
    G == semidefinite(n);
    for i = 1:n
        for j = 1:(i-1)
            if A_local(i,j) == 1
                G(i,i)+G(j,j)-2*G(i,j) <= r^2;
            elseif A_local(i,j) == 0
                G(i,i)+G(j,j)-2*G(i,j) >= r^2 + epsilon;
            end
        end
    end
    sum(sum(G)) == 0;
cvx_end


flag = cvx_optval;

% extracting a solution of correct dimension
[V,L] = eig(G,'vector');
[L, ind] = sort(L,'descend');
V = V(:, ind);
L = diag(L);
X0 = V(:,1:dim)*(L(1:dim, 1:dim))^0.5;

if local == "true"
    % gradient descent
        f = @(x) objective_function(x,A,r,method = method,edges = edges_input, edge_lengths = edge_lengths, backbone = backbone, backbone_length = backbone_length);
        X = fminunc(f,X0);
else
    X = X0;
end
