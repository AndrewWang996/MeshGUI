function [edges, weights, predecessor] = getSpanningTree(meshname)
    
    [V,F] = getMesh(meshname);
    G = meshToGraph(V,F);
    
    anchorIndices = getAnchorIndices(meshname);
    anchorIndex = anchorIndices(1);

    [tree, predecessor] = minspantree(G,'Root',anchorIndex);
    
    edges = orderedDfs(tree, anchorIndex);
    edges = sortrows(edges, 2);
    
    weights = getDistance( V(edges(:,1),:), V(edges(:,2),:) );
    
    if size(predecessor, 1) == 1
        predecessor = transpose(predecessor);
    end
    
end


function dist = getDistance(p1, p2)
    dist = abs(complex(p1(:,1),p1(:,2)) - complex(p2(:,1),p2(:,2)));
end
    