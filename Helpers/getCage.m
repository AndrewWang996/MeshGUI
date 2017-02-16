function cagepts = getCage(meshname)
    cagefiletype = '.mat';
    meshdirectory = 'Meshes/';
    cagefilepath = strcat(meshdirectory, meshname, '/', meshname, cagefiletype);
    if exist(cagefilepath, 'file') == 2
        cagepts = getfield( load(cagefilepath), 'cagepts' );
    else
        cagepts = boundaryPicker(meshfilepath);
        save( cagefilepath, 'cagepts' );
    end
end