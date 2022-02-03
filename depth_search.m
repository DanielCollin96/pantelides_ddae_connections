function connections = depth_search(current,input_node,G,color,connections)
% Depth-first search algorithm to compute all connections in the shifting
% graph. Possibly computes many duplicates.

% Increase length of connection that is currently assembled by adding the
% input node.
current.len = current.len + 1;
current.set(current.len*3-2) = input_node;
color(input_node) = 1;

% Get connected nodes of input node.
out_nodes = find(G.adjacency_matrix(:,input_node));

% Loop over connected uncolored nodes.
for i = 1:length(out_nodes)
    if color(out_nodes(i)) == 0
        
        % Add alternating path to current connection.
        color_temp = color;
        current.set(current.len*3-1) = out_nodes(i);
        color_temp(out_nodes(i)) = 1;
        current.set(current.len*3) = G.assign(out_nodes(i));
        color_temp(out_nodes(i)) = 1;
        
        % Store connection if it is fully assembled.
        if current.len == G.n_assigned
            
            % Increase size of array if necessary.
            if connections.counter > length(connections.array(:,1))
                connections_new = [connections.array; zeros(10000,G.n_assigned*3)];
                connections.array = connections_new;
            end
            connections.array(connections.counter,:) = current.set;
            connections.counter = connections.counter + 1;

        else
            % Recurse from all equation nodes in the current connection.
            connections = depth_search(current,current.set(current.len*3),G,color_temp,connections);
            for j = 1:current.len
                connections = depth_search(current,current.set(j*3-2),G,color_temp,connections);
            end
        end
    end
end

end