function Aeq = getAeq(meshname, vertexIndices, whichKeyframe, anchorIndex)
    V = getMesh(meshname);
    Aeq = zeros( length(vertexIndices) , length(V) );
    for i = 1 : length(vertexIndices)
        Aeq(i,1:length(V)) = getAeqSingle(meshname, vertexIndices(i), whichKeyframe, anchorIndex);
    end
end


function Aeq = getAeqSingle(meshname, vertexIndex, whichKeyframe, anchorIndex)
    V = getMesh(meshname);
    
    keyframe = getKeyframe(meshname, whichKeyframe);
    fz = keyframe.fz;
    
    pathIndices = getPathFromPoint(meshname, anchorIndex, vertexIndex);
    
    Aeq = zeros( 1 , size(V, 1) );
    for i = 2 : length(pathIndices)
        Aeq( pathIndices(i) ) = 1 / conj( fz( pathIndices(i) ) );
    end
end


