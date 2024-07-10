function y = objective_function(x,A,r,options)

% Function represents the objective function in the gradient descent
%
% Input:
% A      Hi-C matrix
% r      radius in which an edge is observed 
% method reconstruction method, options 'unit_ball', 'marble_graph'

arguments
    x
    A
    r
    options.method
    options.edges
    options.edge_lengths
    options.backbone
    options.backbone_length
end

if isfield(options,'method')
    method = options.method;
else
    method = 'unit_ball';
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

n = size(A,1);

if isfield(options,'edge_lengths')
    edge_lengths = options.edge_lengths;
else
    edge_lengths = zeros(n,n);
end

if isfield(options,'edges')
    edges = options.edges;
else
    edges = [];
end

% reset edges for methods which should not use them
if strcmp(method,'marble_graph') || strcmp(method,'unit_ball')
    edges = [];
end

% construct backbone edges
if backbone == "true"
    back_edges = zeros(n-1,2);
    for i = 2:n
        edge_lengths(i-1,i) = backbone_length^2;
        edge_lengths(i,i-1) = backbone_length^2;
        back_edges(i-1,1:2) = [i-1,i];
    end
    edges = [edges;back_edges];
end

A_local=A;

% make sure that edges are not treated as inequalities
for k = 1:size(edges,1)
    i = edges(k,1);
    j = edges(k,2);
    A_local(i,j)=-1;
    A_local(j,i)=-1;
end

y = 0;

for i=1:n
    for j=1:(i-1)
        if A_local(i,j)==1
                y = y + (max(0,norm(x(i,:)-x(j,:))-r))^2;
        elseif A_local(i,j)==0
            y = y + (max(0,-(norm(x(i,:)-x(j,:)))+r))^2;
        else
            y = y + (norm(x(i,:)-x(j,:))-sqrt(edge_lengths(i,j)))^2;
        end
    end
end


