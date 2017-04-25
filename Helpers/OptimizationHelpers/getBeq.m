function Beq = getBeq(meshname, vertexIndices, df_dt)
    Beq = zeros( length(vertexIndices) , 1 );
    for i = 1 : length(vertexIndices)
        Beq(i) = getBeqSingle(meshname, vertexIndices(i), df_dt(i));
    end
end

function BeqSingle = getBeqSingle(meshname, vertexIndex, df_dt)
    [V,F] = getMesh(meshname);
    
    keyframe_1 = getKeyframe(meshname, 1);
    keyframe_2 = getKeyframe(meshname, 2);
    
    pathIndices = getPathFromAnchor(meshname, vertexIndex);
    
    
    fz_1 = keyframe_1.fz;
    fz_2 = keyframe_2.fz;
    logFz_1 = getLogFz(meshname, 1);
    logFz_2 = getLogFz(meshname, 2);
    fzbar_1 = keyframe_1.fzbar;
    
    
    integral_dfz_dt = trapezoidSum( ...
        V, ...
        pathIndices, ...
        @(ind) fz_1(ind) .* ( logFz_2(ind) - logFz_1(ind) ) ...
    );
    
    eta_1 = conj(fz_1) .* fzbar_1;
    
    integral_dfzbar_dt = trapezoidSum( ...
        V, ...
        pathIndices, ...
        @(ind) eta_1(ind) .* conj( (logFz_1(ind) - logFz_2(ind)) ./ fz_1(ind) ) ...
    );
    
    BeqSingle = df_dt - ( integral_dfz_dt + integral_dfzbar_dt );
end


function integralSum = trapezoidSum(vertices, pathIndices, func)
    vert_comp = complex(vertices(:,1), vertices(:,2));
    
    weights = abs(vert_comp(pathIndices(2:end)) - vert_comp(pathIndices(1:end-1)));
    func_avgs = ( func(pathIndices(2:end)) + func(pathIndices(1:end-1)) ) / 2;
    
    integralSum = dot(weights, func_avgs);
end






