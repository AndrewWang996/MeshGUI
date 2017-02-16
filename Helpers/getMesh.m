function [V,F] = getMesh(meshname)
    meshdirectory = 'Meshes/';
    meshfiletype = '.off';

    meshfilepath = strcat(meshdirectory, meshname, '/', meshname, meshfiletype);
    [V,F] = readOff(meshfilepath);

end