function cagepts = getCage(meshname)
    cagefiletype = '.mat';
    meshfiletype = '.off';
    meshdirectory = 'Meshes/';
    cagefilepath = strcat(meshdirectory, meshname, '/', meshname, cagefiletype);
    meshfilepath = strcat(meshdirectory, meshname, '/', meshname, meshfiletype);
    if exist(cagefilepath, 'file') == 2
        cagepts = getfield( load(cagefilepath), 'cagepts' );
    else
        cagepts = boundaryPicker(meshfilepath);
        save( cagefilepath, 'cagepts' );
    end
end