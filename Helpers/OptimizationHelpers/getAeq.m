function Aeq = getAeq(meshname, vertexIndices)
    V = getMesh(meshname);
    Aeq = zeros( length(vertexIndices) , length(V) );
    for i = 1 : length(vertexIndices)
        Aeq(i,1:length(V)) = getAeqSingle(meshname, vertexIndices(i));
    end
    
end


function Aeq = getAeqSingle(meshname, vertexIndex)
    V = getMesh(meshname);
    
    keyframe_0 = getKeyframe(meshname, 1);
    fz_0 = keyframe_0.fz;
    
    pathIndices = getPathFromAnchor(meshname, vertexIndex);
    
    Aeq = zeros( 1 , size(V, 1) );
    for i = 2 : length(pathIndices)
        Aeq(i) = 1 / conj( fz_0( pathIndices(i) ) );
    end
end


function d = dist(z1, z2)
    d = abs(z1 - z2);
end


