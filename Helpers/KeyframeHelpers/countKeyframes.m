function count = countKeyframes(meshname)
    parentDirectory = strcat('Meshes/', meshname, '/');
    keyframeDirectory = strcat(parentDirectory, 'keyframes');
    
    if exist(keyframeDirectory, 'dir') ~= 7
        mkdir(parentDirectory, 'keyframes');
        count = 0;
    else
        D = dir(keyframeDirectory);
        count = length( D(not([D.isdir])) );
    end
end

