function connections = delete_duplicates(connections)
% Deletes all double connections.

if isempty(connections)
    return
end

n_connections = length(connections(:,1));

for i = 1:n_connections
    connections_mod = reshape(connections(i,:),3,[])';
    connections_mod = sortrows(connections_mod);
    connections(i,:) = reshape(connections_mod',1,[]);
end

connections = sortrows(connections);

connections = unique(connections,'rows');

end