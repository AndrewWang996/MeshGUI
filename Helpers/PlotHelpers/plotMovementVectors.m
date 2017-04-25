function handle = plotMovementVectors(vectorsA, vectorsB, scale)
    if nargin < 3
        scale = 0;
    end
    displacement = vectorsB - vectorsA;
    handle = quiver(vectorsA(:,1), vectorsA(:,2), ...
                    displacement(:,1), displacement(:,2), ...
                    scale, ...
                    'color', [0 0 1]);
end