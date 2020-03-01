function comp_num_edges = get_num_edges_components(G)

% This function returns the number of edges in each component of G

G = graph(G);

comps = conncomp(G); % assigns each node in G a component that it belongs to
comp_num_edges = zeros(max(comps),1);

for ii = 1:max(comps) % iterate through each component
    idx = find(comps==ii); % find the nodes that correspond to a single component
    conn_comp = subgraph(G,idx); % isolate the component
    comp_num_edges(ii) = numedges(conn_comp); % count the number of edges in the component
end
