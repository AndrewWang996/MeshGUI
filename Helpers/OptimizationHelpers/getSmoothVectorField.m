function vectorField = getSmoothVectorField(V, F, indices, vectorValues)
    sparseCotMatrix = cotmatrix(V, F);
    f = zeros( size(V,1) );
    A = sparse( size(V,1), size(V,1) );
    b = sparse( size(V,1), 1 );
    options = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex', 'StepTolerance', 1e-15);
    
    vectorField = zeros( size(vectorValues) );
    for dimension = 1:size(vectorValues,1)
        Aeq = sparse(indices, indices, ones(size(indices,1),1), size(V,1));
        beq = vectorValues(dimension, :);
        x = quadprog(sparseCotMatrix, f, A, b, Aeq, beq, options);
        
        vectorField(dimension,:) = x;
    end

end