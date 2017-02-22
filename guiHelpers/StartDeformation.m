function StartDeformation(src,event)
    meshname = getfield( get(gcf, 'UserData'), 'meshname' );
    
    indices = getAnchorIndices(meshname);
    % anchors = 
    display('select start and end for deformation')

end