function accumulatedValues = accumulateAlongEdges(tree, anchorIndex, values)
%% Note that weights also refers to distances between nodes
% tree is a given as a graph without cycles, whose root is anchorIndex
% values is to be accumulated along the edges of the tree

% not really sure what the fastest method is for this...


    T = dfsearch(tree, anchorIndex, 'edgetonew');
    Ta = T(:,1);
    Tb = T(:,2);

    accumulatedValues = zeros( size(values) );
    for i = 1 : size(T,1)
        accumulatedValues( Tb(i) ) = accumulatedValues( Ta(i) ) + values( Tb(i) );
    end
end