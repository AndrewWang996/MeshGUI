function Beq = getBeq(meshname, vertexIndices, df_dt, whichKeyframe, anchorIndex)
    Beq = zeros( length(vertexIndices) , 1 );
    for i = 1 : length(vertexIndices)
        Beq(i) = getBeqSingle(meshname, vertexIndices(i), df_dt(i), whichKeyframe, anchorIndex);
    end
end

function BeqSingle = getBeqSingle(meshname, vertexIndex, df_dt, whichKeyframe, anchorIndex)
    [V,F] = getMesh(meshname);
    
    
    pathIndices = getPathFromPoint(meshname, anchorIndex, vertexIndex);
    
    nKeyframes = countKeyframes(meshname);
    allVertices = zeros(size(V,1), nKeyframes);
    allFz = zeros(size(V,1), nKeyframes);
    allFzbar = zeros(size(V,1), nKeyframes);
    allLogFz = zeros(size(V,1), nKeyframes);

    for i = 1:nKeyframes
        keyframe = getKeyframe(meshname, i);
        allVertices(:,i) = complex(...
            keyframe.Vertices(:,1),...
            keyframe.Vertices(:,2)...
        );
        allFz(:,i) = keyframe.fz;
        allFzbar(:,i) = keyframe.fzbar;
        
        allLogFz(:,i) = getLogFz(meshname, i);
    end
    
    coeffs = getCoefficientsFromControlPoints(allLogFz);
    
    
    function r = func_dfz_dt(ind)
        if whichKeyframe < nKeyframes
            r = coeffs(ind,3) .* allFz(ind, whichKeyframe);
        else
            r = (3 * coeffs(ind,1) + 2 * coeffs(ind,2) + coeffs(ind,3)) ...
                .* allFz(ind, whichKeyframe);
        end
    end

    eta = conj( allFz(:,whichKeyframe) ) .* allFzbar(:,whichKeyframe);
    
    function r = func_dfzbar_dt(ind)
        if whichKeyframe < nKeyframes
            r = eta(ind) .* conj( - coeffs(ind,3) ./ allFz(ind, whichKeyframe) );
        else
            r = eta(ind) .* conj( - (3 * coeffs(ind,1) + 2 * coeffs(ind,2) + coeffs(ind,3) ) ...
                ./ allFz(ind, whichKeyframe) );
        end
    end
    
    integral_dfz_dt = trapezoidSum( ...
        V, ...
        pathIndices, ...
        @func_dfz_dt ...
    );
    
    integral_dfzbar_dt = trapezoidSum( ...
        V, ...
        pathIndices, ...
        @func_dfzbar_dt ...
    );
    
    BeqSingle = df_dt - ( integral_dfz_dt + integral_dfzbar_dt );
end


function integralSum = trapezoidSum(vertices, pathIndices, func)
    complex_vertices = complex(vertices(:,1), vertices(:,2));
    
    weights = abs(complex_vertices(pathIndices(2:end)) - complex_vertices(pathIndices(1:end-1)));
    func_avgs = ( func(pathIndices(2:end)) + func(pathIndices(1:end-1)) ) / 2;
    
    integralSum = dot(weights, func_avgs);
end






