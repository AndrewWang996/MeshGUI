function keyframe = getKeyframe(meshname, whichKeyframe)
    keyframeFile = strcat('Meshes/', meshname, '/keyframes/keyframe_', int2str(whichKeyframe));
    if exist(keyframeFile, 'file') > 0
        savedStruct = load(keyframeFile);
        keyframe = savedStruct.keyframe;
    else
        error('keyframe %i for mesh %s not found', whichKeyframe, meshname);
    end

end