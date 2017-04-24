function [endNodes, weights, predecessor] = getSpanningTree(meshname)
%%% TODO: Enforce user selection of anchor.
% Currently, MATLAB's minspantree sets node 1 as the root
% This may not be optimal, so in the future we want the 1st anchor
% as root.
%
% Right now, this is a little annoying to do, and not strictly necessary. 
    
    [V,F] = getMesh(meshname);
    G = meshToGraph(V,F);

    [tree, predecessor] = minspantree(G);
    endNodes = tree.Edges.EndNodes;
    weights = tree.Edges.Weight;
    
    if size(predecessor, 1) == 1
        predecessor = transpose(predecessor);
    end
    
end

function dist = getDistance(p1, p2)
    dist = abs(p1 - p2);
end

    
    