function [endNodes, weights, predecessor] = getSpanningTree(meshname)
    
    [V,F] = getMesh(meshname);
    G = meshToGraph(V,F);
    
    anchorIndices = getAnchorIndices(meshname);
    anchorIndex = anchorIndices(1);

    [tree, predecessor] = minspantree(G,'Root',anchorIndex);
    endNodes = tree.Edges.EndNodes;
    weights = tree.Edges.Weight;
    
    if size(predecessor, 1) == 1
        predecessor = transpose(predecessor);
    end
    
end


    
    