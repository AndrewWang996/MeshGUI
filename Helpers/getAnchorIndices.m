function indices = getAnchorIndices(meshname)
    
    anchorfiletype = '.mat';
    meshfiletype = '.off';
    meshdirectory = 'Meshes/';
    meshfilepath = strcat(meshdirectory, meshname, '/', meshname, meshfiletype);
    anchorfilepath = strcat(meshdirectory, meshname, '/', 'anchor', anchorfiletype);
    
    if exist(anchorfilepath, 'file') == 2
        indices = getfield( load(anchorfilepath), 'indices' );
    else
        indices = anchorPicker(meshfilepath);
        save( anchorfilepath, 'indices' );
    end

end


