function saveKeyframe(meshname, keyframe)
%%% Saves a struct keyframe consisting of
% - Vertices
% - Faces
% - fz
% - fzbar
%%% into a folder

    parentDirectory = strcat('Meshes/', meshname, '/');
    keyframeDirectory = strcat(parentDirectory, 'keyframes');
    
    if exist(keyframeDirectory, 'dir') ~= 7
        mkdir(parentDirectory, 'keyframes');
    end
    
    numKeyframes = countKeyframes(meshname);
    keyframeFile = strcat(keyframeDirectory, '/keyframe_', int2str(numKeyframes + 1));
    save(keyframeFile, 'keyframe');

end

