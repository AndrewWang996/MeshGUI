function save_deta_dt(meshname, whichKeyframe, deta_dt)
    
    parentDirectory = strcat('Meshes/', meshname, '/');
    keyframeDirectory = strcat(parentDirectory, 'deta_dt');
    
    if exist(keyframeDirectory, 'dir') ~= 7
        mkdir(parentDirectory, 'deta_dt');
    end
    
    deta_dt_file = strcat(keyframeDirectory, '/deta_dt_', int2str(whichKeyframe));
    save(deta_dt_file, 'deta_dt');

end