function Beq = getBeq(meshname, vertexIndices, df_dt, whichKeyframe, anchorIndex)
    Beq = zeros( length(vertexIndices) , 1 );
    for i = 1 : length(vertexIndices)
        Beq(i) = getBeqSingle(meshname, vertexIndices(i), df_dt(i), whichKeyframe, anchorIndex);
    end
end

function BeqSingle = getBeqSingle(meshname, vertexIndex, df_dt, whichKeyframe, anchorIndex)
    [V,F] = getMesh(meshname);
    
    if whichKeyframe == 1
        A = 1;
        B = 2;
    else
        A = whichKeyframe - 1;
        B = whichKeyframe;
    end
    
    pathIndices = getPathFromPoint(meshname, anchorIndex, vertexIndex);
    
    keyframe_A = getKeyframe(meshname, A);
    keyframe_B = getKeyframe(meshname, B);
    
    
    
    fz_A = keyframe_A.fz;
    fz_B = keyframe_B.fz;
    logFz_A = getLogFz(meshname, A);
    logFz_B = getLogFz(meshname, B);
    fzbar_A = keyframe_A.fzbar;
    fzbar_B = keyframe_B.fzbar;
    
    
    function r = func_dfz_dt(ind)
        if whichKeyframe == 1
            r = fz_A(ind) .* (logFz_B(ind) - logFz_A(ind));
        else
            r = fz_B(ind) .* (logFz_B(ind) - logFz_A(ind));
        end
    end

    eta = conj(fz_A) .* fzbar_A;
    
    function r = func_dfzbar_dt(ind)
        if whichKeyframe == 1
            r = eta(ind) .* conj( (logFz_A(ind) - logFz_B(ind)) ./ fz_A(ind) );
        else
            r = eta(ind) .* conj( (logFz_A(ind) - logFz_B(ind)) ./ fz_B(ind) );
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






