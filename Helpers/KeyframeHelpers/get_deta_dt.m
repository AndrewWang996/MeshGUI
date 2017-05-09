function deta_dt = get_deta_dt(meshname, whichKeyframe)
    deta_dt_file = strcat('Meshes/', meshname, '/deta_dt/deta_dt_', int2str(whichKeyframe), '.mat');
    if exist(deta_dt_file, 'file') > 0
        savedStruct = load(deta_dt_file);
        deta_dt = savedStruct.deta_dt;
    else
        V = getMesh(meshname);
        deta_dt = ones( size(V,1), 1 );
        fprintf('deta_dt not found for keyframe %i, defaulting to zero vector.\n', whichKeyframe);
    end
end