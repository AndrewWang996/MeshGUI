function Beq = getBeq(meshname, vertexIndices, df_dt)
    Beq = zeros( length(vertexIndices) , 1 );
    for i = 1 : length(vertexIndices)
        Beq(i) = getBeqSingle(meshname, vertexIndices(i), df_dt(i));
    end
end

function BeqSingle = getBeqSingle(meshname, vertexIndex, df_dt)

    cage = getCage(meshname);
    [V,F] = getMesh(meshname);
    
    keyframe_0 = getKeyframe(meshname, 1);
    
    phi_0 = keyframe_0.phi;
    
    pathIndices = getPathFromAnchor(meshname, vertexIndex);
    pathVertices = V(pathIndices,:);
    
    integral_dfz_dt = integral(...
        @(z) fz(cage, z, phi_0) .* (logFz(meshname, 2, z) - logFz(meshname, 1, z)), ...
        complex(pathVertices(1,1), pathVertices(1,2)), ...
        complex(pathVertices(end,1), pathVertices(end,2)), ...
        'Waypoints', complex(pathVertices(2:end-1,1), pathVertices(2:end-1,2)) ...
    );

    integral_dfzbar_dt = integral(...
         @(z) eta(cage, z, keyframe_0.phi, keyframe_0.psi) ...
            .* conj( (logFz(meshname, 1, z) - logFz(meshname, 2, z)) ./ fz(cage, z, phi_0) ), ...
        complex(pathVertices(1,1), pathVertices(1,2)), ...
        complex(pathVertices(end,1), pathVertices(end,2)), ...
        'Waypoints', complex(pathVertices(2:end-1,1), pathVertices(2:end-1,2)) ...
    );
    
    BeqSingle = df_dt - ( integral_dfz_dt + integral_dfzbar_dt );
end


function r = logFz(meshname, whichKeyframe, z)
    [V,F] = getMesh(meshname);
    cage = getCage(meshname);
    keyframe = getKeyframe(meshname, whichKeyframe);
    phi = keyframe.phi;
    
    vertexIndices = getIndex(reshape(z, length(z), 1), V(:,1:2));
    
    r = anchorLogFz(meshname, 1) * ones( length(z), 1 );
    
    for i = 1:length(z)
        vertexIndex = vertexIndices(i);
        pathIndices = getPathFromAnchor(meshname, vertexIndex);
        pathVertices = complex( V(pathIndices,1), V(pathIndices,2) );

        errorTerm = log( fz(cage, z(i), phi) ./ fz(cage, complex(V(vertexIndex,1), V(vertexIndex,2)), phi) );
        if length(pathIndices) < 2
            r(i) = r(i) + errorTerm;
            continue
        end
        fracList = fz(cage, pathVertices(2:end), phi) ./ fz(cage, pathVertices(1:end-1), phi);
        summation = sum( log( fracList ) );
        
        r(i) = r(i) + summation + errorTerm;
    end
    
    r = reshape(r, 1, length(r));
end

function r = anchorLogFz(meshname, whichKeyframe)
    anchorIndices = getAnchorIndices(meshname);
    anchorIndex = anchorIndices(1);
    keyframe_0 = getKeyframe(meshname, 1);
    keyframe = getKeyframe(meshname, whichKeyframe);
    fz_0 = keyframe_0.fz;
    fz = keyframe.fz;
    r = log(fz(anchorIndex) / fz_0(anchorIndex)) + log( fz_0(anchorIndex) );
end

function r = fz(cage, z, phi)
    D = derivativesOfCauchyCoord(cage, z);
    r = D*phi;
    r = reshape(r, 1, length(r));
end

function r = eta(cage, z, phi, psi)
    D = derivativesOfCauchyCoord(cage, z);
    fz = D * phi;
    fzbar = conj( D * psi );
    
    r = conj(fz) .* fzbar;
    r = reshape(r, 1, length(r));
end






