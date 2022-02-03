% This code provides a comparison between a naive and simple depth-first
% search method and the new method based on the enumeration of the spanning
% trees in the connection graph.
% Author: Daniel Collin, Technische Universität Berlin,
% daniel.collin@web.de.

% Specify system size and the structure of the shifting graph, i.e.,
% equ_var_structure(i,j) = 1 if v_j appears in F_i.

n = 4;
% equ_var_structure = triu(ones(n,n),0);
% equ_var_structure = equ_var_structure(:,1:n-1);
% equ_var_structure(n,:) = ones(1,n-1);

% equ_var_structure = diag(ones(n,1)) + diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
% equ_var_structure = equ_var_structure(:,1:n-1);
% equ_var_structure(n,:) = ones(1,n-1);

equ_var_structure = ones(n,n-1);


equ_var_structure




% Create shifting graph and matching stored in the vector assign, i.e.,
% assign(i) = j if node with index i is matched to equation F_i.
shifting_graph = struct;
shifting_graph.adjacency_matrix = [zeros(n,n), equ_var_structure; equ_var_structure', zeros(n-1,n-1)];                             
shifting_graph.assign = [zeros(1,n), 1:n-1]';
shifting_graph.n_assigned = n-1;

% Here, it is assumed for simplicity that the last node F_n is exposed.
exposed_node = n;

%% Depth-first search algorithm
tic

% Vector that colors visited nodes.
color = zeros(2*n-1,1);

% Structure to store connections.
connections = struct;
connections.array = zeros(10000,3*(n-1));
connections.counter = 1;

% Structure to store connection that is currently assembled.
current = struct;
current.set = zeros(3*(n-1),1);
current.len = 0;

% Start recursion.
connections = depth_search(current,exposed_node,shifting_graph,color,connections);

% Set of all connections with duplicates.
P_with_duplicates = connections.array(1:connections.counter-1,:);

% Delete duplicates
P_dfs = delete_duplicates(P_with_duplicates);
toc

fprintf('Number of connections computed by DFS before deleting duplicates: %i\n',size(P_with_duplicates,1));
fprintf('Number of connections computed by DFS after deleting duplicates: %i\n',size(P_dfs,1));



%% New method based on the enumeration of the spanning trees in the connection graph
tic

% Create structures for enumeration algorithm.
H = struct;
T = struct;
dl = struct;
L = struct;


% Construct connection graph.
H.adj = zeros(n,n);
H.n = n;

% Loop over or all variable nodes.
for i = n+1:2*n-1
    
    % Get connected nodes.
    out_nodes = find(shifting_graph.adjacency_matrix(:,i));
    
    % Loop over connected nodes.
    for j = 1:length(out_nodes)
        
        % Add alternating paths and variables.
        if out_nodes(j) ~= shifting_graph.assign(i)
            H.adj(out_nodes(j),shifting_graph.assign(i)) = i;
        end
    end
end

H.sp_trees = zeros(1000,3*(H.n-1));
H.counter = 1;


% Create nodes for a doubly linked list and store them in dl, and create
% an array dl_pointer where each entry points to the doubly linked node of
% the corresponding edge.
nEdges = nnz(H.adj);
dl.nodes = cell(nEdges,1);
[ii,jj] = find(H.adj);
dd = 1:nEdges;
H.dl_pointer = sparse(ii,jj,dd,H.n,H.n);

for k = 1:nEdges
    dl.nodes{k} = dlnode([ii(k),jj(k)]);
end

% Initial tree contains only the root.
T.adj = spalloc(H.n,H.n,H.n-1);
T.indicator = zeros(H.n,1);
T.indicator(exposed_node) = 1;

% Initialize vectors for preorder numbers.
L.P = zeros(H.n,1);
L.P(exposed_node) = 1;

L.H = zeros(H.n,1);
L.preorder_counter = 1;

% Create initial F as doubly linked list of outgoing edges from T. F is 
% used as a stack. dl.F always points to the last added element.
node_id = find(H.adj(exposed_node,:));
dl.F = dl.nodes{H.dl_pointer(exposed_node,node_id(1))};

for i = 2:length(node_id)
    insertAfter(dl.nodes{H.dl_pointer(exposed_node,node_id(i))},dl.F);
    dl.F = dl.F.Next;
end

% Compute spanning trees.
[H,~] = grow(H,T,dl,L);


% Shorten array to the needed size.
P_grow = H.sp_trees(1:H.counter-1,:);
toc

fprintf('Number of connections computed by GROW: %i\n',size(P_grow,1));


