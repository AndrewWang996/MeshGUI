function edges = orderedDfs(tree, root)
    numVertices = numnodes(tree);
    edges = dfsearch(tree, root, 'edgetonew');
    
    visited = zeros(numVertices);
    visited(root) = 1;
    
    for i = 1:length(edges)
        v1 = edges(i,1);
        v2 = edges(i,2);
        if ~ visited(v1)
            edges(i,:) = [v2 v1];
            visited(v1) = 1;
        elseif ~ visited(v2)
            edges(i,:) = [v1 v2];
            visited(v2) = 1;
        else
            disp('something is wrong with your dfs');
        end
    end

end