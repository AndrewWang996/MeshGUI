function deleted = deleteKeyframe(meshname, whichKeyframe)
    fileDirectory = strcat('Meshes/', meshname, '/keyframes/');
    keyframeFile = strcat(fileDirectory, 'keyframe_', int2str(whichKeyframe), '.mat');
    display(keyframeFile)
    if exist(keyframeFile, 'file') > 0
        delete(keyframeFile);
        for i = whichKeyframe + 1 : countKeyframes(meshname)
            nextKeyframeFile = strcat(fileDirectory, 'keyframe_', int2str(i), '.mat');
            if exist(nextKeyframeFile, 'file') <= 0
                break
            end
            movefile(nextKeyframeFile, keyframeFile);
            keyframeFile = nextKeyframeFile;
        end
        deleted = true;
    else
        deleted = false;
    end

end