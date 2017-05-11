function s = getInformation(meshname)
% return a struct s with information:
% - anchorIndex ...
% - vertices ...
% - faces ...
% - allVertices ...
% - allFz ... 
% - allFzBar ... 
% - allLogFz ...
% - allEta ...
% - all_deta_dt ...
% - endNodes ...

    s = struct;
    
    numKeyframes = countKeyframes(meshname);
    [vertices, faces] = getMesh(meshname);
    s.vertices = vertices;
    s.faces = faces;
    
    numVertices = size(vertices,1);
    
    allVertices = zeros(numVertices, numKeyframes);
    allFz = zeros(numVertices, numKeyframes);
    allFzBar = zeros(numVertices, numKeyframes);
    all_deta_dt = zeros(numVertices, numKeyframes);
    
    
    

    for whichKeyframe = 1:numKeyframes
        keyframe = getKeyframe(meshname, whichKeyframe);
        allVertices(:,whichKeyframe) = complex(...
            keyframe.Vertices(:,1),...
            keyframe.Vertices(:,2)...
        );
        allFz(:,whichKeyframe) = keyframe.fz;
        allFzBar(:,whichKeyframe) = keyframe.fzbar;

        deta_dt = get_deta_dt(meshname, whichKeyframe);
        all_deta_dt(:,whichKeyframe) = deta_dt;
    end

    allEta = conj(allFz) .* allFzBar;
    
    s.allVertices = allVertices;
    s.allFz = allFz;
    s.allFzBar = allFzBar;
    s.allEta = allEta;
    s.all_deta_dt = all_deta_dt;
    
    anchorIndices = getAnchorIndices(meshname);
    anchorIndex = anchorIndices(1); % only use the first anchor
    s.anchorIndex = anchorIndex;
    
    edges = getSpanningTree(meshname);
    endNodes = [edges(1:anchorIndex-1,:);...
        anchorIndex, anchorIndex;...
        edges(anchorIndex:end,:)];
    
    s.endNodes = endNodes;


    angleDiffs = accumulateAlongEdges(...
        graph(endNodes(:,1), endNodes(:,2)),...
        anchorIndex,...
        angle( allFz(endNodes(:,2),:) ./ allFz(endNodes(:,1),:) ),...
        angle( allFz(anchorIndex,:) )...
    );

    allLogFz = complex( log(abs(allFz)), angleDiffs );    
    s.allLogFz = allLogFz;

end