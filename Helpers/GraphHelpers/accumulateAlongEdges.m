function accumulatedValues = accumulateAlongEdges(tree, anchorIndex, values, initial_values)
%% Note that weights also refers to distances between nodes
% tree is a given as a graph without cycles, whose root is anchorIndex
% values is to be accumulated along the edges of the tree

% not really sure what the fastest method is for this...


    T = orderedDfs(tree, anchorIndex);

    accumulatedValues = zeros( size(values) );
    
    if nargin >= 4
        accumulatedValues(anchorIndex,:) = initial_values;
    end
    for i = 1:length(T)
        accumulatedValues(T(i,2),:) = accumulatedValues(T(i,1),:) + values(T(i,2),:);
    end
end