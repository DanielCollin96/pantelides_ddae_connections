function [G,L] = grow(G,T,dl,L)
% This is an implementation of the algorithm presented in the paper
% "Finding all spanning trees of directed and undirected graphs" from
% Harold N. Gabow and Eugene W. Myers (1978).
% grow computes all spanning trees of G rooted at r containing the tree T.
% For a detailed explanation of the algorithm we refer
% to the paper.
%
% Input:
% - G: directed graph G with the following properties:
%   - G.adj: adjacency matrix with g_ij being the edge weight (or the
%     number k of the variable x_k of an alternating path (F_i,x_k,F_j)).
%   - G.n: number of nodes of G.
%   - G.dl_pointer: array that points on the position of the doubly linked 
%     node of each edge in dl.nodes such that 
%     dl_pointer(i,j) = k if dl{k} is the doubly linked node of the edge 
%     (v_i,v_j).
%   - G.sp_trees: array that stores computed spanning trees row-wise.
%   - G.counter: number of computed spanning trees - 1.
% - T: subgraph of G that represents the spanning tree that is being built
%   with the following properties:
%   - T.adj: adjacency matrix with t_ij being the edge weight (or the
%     number k of the variable x_k of an alternating path (F_i,x_k,F_j)).
%   - T.indicator: vector that indicates which nodes of G are part of T.
% - dl: structure of doubly linked list nodes with the following
%   properties:
%   - dl.nodes: cell array that stores the doubly linked nodes belonging to
%     the edges of the graph G.
%   - dl.F: node that was most recently added to the doubly linked list F.
%     F is managed as a stack and contains all edges that can be addded to
%     T.
% - L: directed graph that represents the spanning tree that has been
%   computed last with the following properties:
%   - L.adj: adjacency matrix with l_ij being the edge weight (or the
%     number k of the variable x_k of an alternating path (F_i,x_k,F_j)).
%   - L.P: vector that contains preorder numbers of the vertices of L,
%     neccessary for the bridge test. Note that for a better efficiency,
%     preorder labeling of trees is done as the trees are grown and does
%     not correspond to the exact preorder numbers as defined in
%     literature. More specifically, several vertices can have the same
%     preorder numbers. However, this does not affect the functionality of
%     the bridge test.
%   - L.H: vector that contains the highest preorder number of the 
%     descendents of the vertices of L, i.e. H(v) is the highest preorder 
%     number of a descendent of v in L.
%   - L.preorder_counter: integer that stores current number of
%     preordering.


nEdges = nnz(T.adj);

%% Save spanning tree

if nEdges == G.n - 1
    L.adj = T.adj;
    
    % Increase array size if necessary
    if G.counter > size(G.sp_trees,1)
        G.sp_trees = [G.sp_trees; zeros(1000,3*(G.n-1))];
    end
    
    % Save edges of T in G.sp_trees
    [E1,E2,var] = find(T.adj);
    E = [E1,var,E2]';
    G.sp_trees(G.counter,:) = reshape(E,1,3*(G.n-1));
    G.counter = G.counter + 1;
    return
end



%% If spanning tree is not complete, iteratively add new edges until a bridge is processed

% Boolean that indicates if edge e is a bridge 
b = false;

% Increase preorder number 
L.preorder_counter = L.preorder_counter + 1;

% Create FF, a local list of edges that have been processed and deleted already
FF = dlnode.empty;
FF_weights = zeros(nnz(G.adj),1);
sizeFF = 0;

while b == false % Until bridge is found
    %% New tree edge
    
    % Pop an edge e=(u,v) from F
    e = dl.F.Data;
    e_node = dl.F;
    dl.F = dl.F.Prev;
    removeNode(e_node);
    
    u = e(1);
    v = e(2);
    
    % Add e to T
    T.adj(u,v) = G.adj(u,v);
    T.indicator(v) = 1;
    
    % Label vertex with current preorder number
    L.P(v) = L.preorder_counter;
    

    %% Update F
    
    % Get all ingoing edges to v except e=(u,v)
    in_node_id = find(G.adj(:,v));
    in_node_id = in_node_id(in_node_id ~= u);
    
    % Remove each edge (w,v), w\in T, from F
    for i = 1:length(in_node_id)
        if T.indicator(in_node_id(i)) == 1
            % Remove edge from F but do not destroy the values of its
            % links to be able to restore it at the same place later
            prevNode = dl.nodes{G.dl_pointer(in_node_id(i),v)}.Prev;
            nextNode = dl.nodes{G.dl_pointer(in_node_id(i),v)}.Next;
            % Update links of previous and next node
            if ~isempty(prevNode)
                prevNode.Next = nextNode;
            end
            if ~isempty(nextNode)
                nextNode.Prev = prevNode;
            else
                % if the current node is the one that was added last to the
                % stack, update the pointer to F
                if isequal(dl.nodes{G.dl_pointer(in_node_id(i),v)},dl.F)
                    dl.F = prevNode;
                end
            end
        end
    end
    
    % Get all outgoing edges from v
    out_node_id = find(G.adj(v,:));
    
    % Push each edge (v,w), w\notin T, onto F
    for i = 1:length(out_node_id)
        if T.indicator(out_node_id(i)) == 0
            if ~isempty(dl.F)
                insertAfter(dl.nodes{G.dl_pointer(v,out_node_id(i))},dl.F);
            end
            % Update last added node of F
            dl.F = dl.nodes{G.dl_pointer(v,out_node_id(i))};
        end
    end
    

    
    %% Recurse
    [G,L] = grow(G,T,dl,L);
    
    % Update H with the highest preorder number of a descendent of v
    node_id = find(L.adj(v,:));
    
    if ~isempty(node_id)
        % Create array to store the number of descendents of v for each
        % outgoing edge
        descendents = zeros(length(node_id),1);
        for i = 1:length(node_id)
            descendents(i) = L.H(node_id(i));
        end
        L.H(v) = max(descendents);
    else
        % if there is no outgoing edge, the highest preorder number of a
        % descendent is the preorder number of the vertex itself
        L.H(v) = L.P(v);
    end
        
    
    
    %% Restore F
    
    % Pop each edge (v,w), w\notin T, from F
    for i = length(out_node_id):-1:1
        if T.indicator(out_node_id(i)) == 0
            % Update pointer to F
            if isequal(dl.nodes{G.dl_pointer(v,out_node_id(i))},dl.F)
                dl.F = dl.F.Prev;
            end
            removeNode(dl.nodes{G.dl_pointer(v,out_node_id(i))});
        end
    end
    
    % Restore each edge (w,v), w\in T, in F
    for i = length(in_node_id):-1:1
        if T.indicator(in_node_id(i)) == 1
            % Insert at the place indicated by the values of the links
            prevNode = dl.nodes{G.dl_pointer(in_node_id(i),v)}.Prev;
            nextNode = dl.nodes{G.dl_pointer(in_node_id(i),v)}.Next;
            if ~isempty(prevNode)
                prevNode.Next = dl.nodes{G.dl_pointer(in_node_id(i),v)};
            end

            if ~isempty(nextNode)
                nextNode.Prev = dl.nodes{G.dl_pointer(in_node_id(i),v)};
            else
                % Update pointer to F if node is inserted at the end of the
                % stack
                if isequal(prevNode,dl.F)
                    dl.F = dl.nodes{G.dl_pointer(in_node_id(i),v)};
                end
            end
        end
    end  
    
    %% Delete e
    
    % Delete e=(u,v) from T
    
    T.adj(u,v) = 0;
    T.indicator(v) = 0;
    
    % Add e to FF
    
    if ~isempty(FF)
        insertAfter(dl.nodes{G.dl_pointer(u,v)},FF);
        FF = dl.nodes{G.dl_pointer(u,v)};
    else
        FF = dl.nodes{G.dl_pointer(u,v)};
    end
    
    sizeFF = sizeFF + 1;
    
    % Delete e from G but preserve its weight to be able to restore it
    % later
    FF_weights(sizeFF) = G.adj(u,v);
    G.adj(u,v) = 0;
    
    
    %% Bridge test
    % Preordering of L stored in P and H was already computed.
    % w is a descendent of v if and only if P(v) \leq P(w) \leq H(v).
    % If there is an edge (w,v), where w is not a descendent of v in L,
    % then b <- false, else b <- true.
    
    % Get all edges that are ingoing to v in G
    in_node_id = find(G.adj(:,v));
    b = true;
    for i = 1:length(in_node_id)
        % If an edge (w,v) is found in G, s.t. (u,v) is not a bridge, the
        % loop can be stopped
        if L.P(v) > L.P(in_node_id(i)) || L.P(in_node_id(i)) > L.H(v)
            b = false;
            break
        end
    end
end % while b == false



%% Reconstruct G

% Create arrays ii and jj for the indices i and j, respectively, for the 
% entries g_ij of the adjacency matrix of G that have to be reconstructed
ii = zeros(sizeFF,1);
jj = zeros(sizeFF,1);

% Shorten the vector of weights of the edges in FF. These will be the
% values of the entries of the adjacency matrix.
FF_weights = FF_weights(1:sizeFF);

while ~isempty(FF)
    % Pop an edge from FF
    e_node = FF;
    FF = FF.Prev;
    removeNode(e_node);
    
    % Push this edge onto F
    if ~isempty(dl.F)
        insertAfter(e_node,dl.F);
        dl.F = dl.F.Next;
    else
        dl.F = e_node;
    end
    
    % Store vertices of edge as row and column index for the adjacency matrix
    e = e_node.Data;
    ii(sizeFF) = e(1);
    jj(sizeFF) = e(2);
    sizeFF = sizeFF - 1;
end

% Add all collected indices as entries to the adjacency matrix of G
G.adj = G.adj + sparse(ii,jj,FF_weights,G.n,G.n);

end