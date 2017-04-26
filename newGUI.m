
addpath Helpers;
addpath Helpers/KeyframeHelpers;
addpath Helpers/PlotHelpers;
addpath Helpers/GUIHelpers;
addpath Helpers/OptimizationHelpers;
addpath Helpers/GraphHelpers;

addpath WeightFunctions;

addpath Meshes;

axis equal;



meshname = 'vert_bar';
[V,F] = getMesh(meshname);
cagepts = getCage(meshname);

set(gcf, 'UserData', struct('meshname', meshname));
set(gcf, 'Position', [0 200 800 500]);

simple_deform(V, F, meshname)




function simple_deform(varargin)

    V = varargin{1};
    % face indices of mesh
    F = varargin{2};

    meshname = varargin{3};
    % control vertices of skeleton
    C = V;
    % what indices of V do C correspond to
    I = 1:length(V);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set default parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Be sure that control vertices are in 2D
    if(size(C,2) >= 3)
        C = C(:,1:2);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Prepare output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global g_Deform;
    gid = numel(g_Deform)+1;
    % store indices I
    g_Deform(gid).indices = I;
    % keep track of control positions at mouse down
    g_Deform(gid).new_C = [];
    g_Deform(gid).update_positions = @update_positions;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Set up figure with enough subplots for any additional visualizations
    % clear current figure
    clf

    % plot the original mesh
    g_Deform(gid).tsh = trisurf(F,V(:,1),V(:,2),zeros(size(V,1),1), ...
        'FaceColor','interp');
    % 2D view
    view(2);
    axis equal
    axis manual
    hold on;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cache fz, fzbar for the identity function f(z) = z
    g_Deform(gid).fz = ones( length(V), 1 );
    g_Deform(gid).fzbar = zeros( length(V), 1 );
    g_Deform(gid).phi = complex(V(:,1), V(:,2));
    g_Deform(gid).psi = zeros( length(V), 1 );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveKeyframeButton = uicontrol(gcf,'Style','pushbutton',...
    'String','Save Keyframe',...
    'Position',[50 0 90 20],...
    'Callback', @SaveKeyframe);

showKeyframeButton = uicontrol(gcf,'Style','slider',...
    'String','Save Keyframe',...
    'Position',[500 0 120 20]);
addlistener(showKeyframeButton, ...
    'ContinuousValueChange', @ShowKeyframe);

showAnimation = uicontrol(gcf,'Style','pushbutton',...
    'String','Show Animation',...
    'Position',[200 0 90 20],...
    'Callback',@ShowAnimation);

startDeformation = uicontrol(gcf,'Style','pushbutton',...
    'String','Deform',...
    'Position',[300 0 90 20],...
    'Callback',@StartDeformation);

setVelocity = uicontrol(gcf,'Style','pushbutton',...
    'String','Set Velocity',...
    'Position',[400 0 90 20],...
    'Callback',@SetVelocity);

    function StartDeformation(src,event)
        disp('click pairs of points, 1st on the shape, 2nd on the desired new location')

        [vertices, faces] = getMesh(meshname);
        cagePts = getCage(meshname);

        figure
        trimesh(faces, vertices(:,1), vertices(:,2), 'color', 'k');
        [ptsX, ptsY] = getpts;
        close(gcf)

        complexPts = complex(ptsX, ptsY);

        ptsFrom = complexPts(1:2:length(complexPts), :);
        ptsTo = complexPts(2:2:length(complexPts), :);

        indices = getIndex(ptsFrom, vertices);
        anchorIndices = getAnchorIndices(meshname);
        anchorPositions = vertices(anchorIndices, 1) + 1i * vertices(anchorIndices, 2);

        [newVerticesComplex, fz, fzbar, phi, psi] = deformBoundedDistortion([indices; anchorIndices], [ptsTo; anchorPositions], vertices, faces, cagePts);
        newVertices = [real(newVerticesComplex), imag(newVerticesComplex)];

        set(g_Deform(gid).tsh, 'Vertices', newVertices);
        g_Deform(gid).fz = fz;
        g_Deform(gid).fzbar = fzbar;
        g_Deform(gid).phi = phi;
        g_Deform(gid).psi = psi;
    end

    hold off;
    
    return
    
    
    function SetVelocity(src,event)
        disp( strcat('Select pairs of points. First point a vertex on mesh.', ...
            ' Second point is endpoint of velocity vector.') );

        [vertices, faces] = getMesh(meshname);

        figure
        trimesh(faces, vertices(:,1), vertices(:,2), 'color', 'k');
        axis equal
        axis manual
        while 1
            [ptsX, ptsY] = getpts; 
            if mod(size(ptsY, 1), 2) ~= 0
                disp( strcat('please remember to select an ', ...
                    'even number of anchor points.\n') );
                continue;
            end
            break;
        end
        close(gcf)

        complexPts = complex(ptsX, ptsY);

        ptsFrom = complexPts(1:2:length(complexPts), :);
        ptsTo = complexPts(2:2:length(complexPts), :);

        indices = getIndex(ptsFrom, vertices);
        hold on;
        plotMovementVectors( ...
            vertices(indices,1:2), ...
            horzcat(real(ptsTo), imag(ptsTo)), ...
            1 ...
        );
        
        
        velocities = ptsTo - complex(vertices(indices,1), vertices(indices,2));
        [H,Aeq,beq] = poseOptimizationProblem(meshname, indices, velocities);
        deta_dt = solveOptimizationProblem(H,Aeq,beq);
        
        
        plotMovementVectors( ...
            vertices, ...
            vertices + horzcat( real(deta_dt), imag(deta_dt) ) ...
        );
        hold off;
    end

    function SaveKeyframe(src, event)
        copyOfCurrent = {};
        copyOfCurrent.Vertices = g_Deform(gid).tsh.Vertices;
        copyOfCurrent.Faces = g_Deform(gid).tsh.Faces;
        copyOfCurrent.fz = g_Deform(gid).fz;
        copyOfCurrent.fzbar = g_Deform(gid).fzbar;
        copyOfCurrent.phi = g_Deform(gid).phi;
        copyOfCurrent.psi = g_Deform(gid).psi;
        
        saveKeyframe(meshname, copyOfCurrent);
        fprintf('Saved 1 new keyframe. %d total keyframes.\n', countKeyframes(meshname));
    end

    function ShowKeyframe(src,event)
        numKeyframes = countKeyframes(meshname);
        if numKeyframes == 0
            return 
        end
        position = src.Value;
        whichKeyframe = round( (numKeyframes - 1) * position + 1);
        
        keyframe = getKeyframe(meshname, whichKeyframe);
        set(g_Deform(gid).tsh,'Vertices',keyframe.Vertices(:,1:2));
    end
    

    function ShowAnimation(src,event)
        
        % for debugging purposes
        showMovementVectors = false;
        
        numKeyframes = countKeyframes(meshname);
        allVertices = zeros(size(C,1), numKeyframes);
        allFz = zeros(size(C,1), numKeyframes);
        allFzbar = zeros(size(C,1), numKeyframes);
        [vertices, faces] = getMesh(meshname);
        
        for whichKeyframe = 1:numKeyframes
            keyframe = getKeyframe(meshname, whichKeyframe);
            allVertices(:,whichKeyframe) = complex(...
                keyframe.Vertices(:,1),...
                keyframe.Vertices(:,2)...
            );
        
            allFz(:,whichKeyframe) = keyframe.fz;
            allFzbar(:,whichKeyframe) = keyframe.fzbar;
        end
        
        % 0) set up spanning tree of mesh
        % + other precomputation        
        [endNodes, distances, predecessor] = getSpanningTree(meshname);
        anchorIndices = getAnchorIndices(meshname);
        anchorIndex = anchorIndices(1); % only use the first anchor
        edgeVectors = complex(...
            vertices(endNodes(:,1),1) - vertices(endNodes(:,2),1),...
            vertices(endNodes(:,1),2) - vertices(endNodes(:,2),2)...
        );
        dfsEdges = dfsearch( graph(endNodes(:,1), endNodes(:,2)), anchorIndex , 'edgetonew' );
        
            
        % 1) do some magic log stuff
        g = complex( log(abs(allFz)), angle(allFz(end)) + cumsum(angle(allFz./circshift(allFz, 1))) );
            
        
        
        function interpF = interpolate(wt, weightFunc)
            weight = weightFunc(numKeyframes, wt);
            
            % 2) interpolate fz, eta, fzbar
            % TODO: define fWtFun(wt) as linear matrix
            % This is already given in BDH code...
            interpFz = exp( g * weight );         
            interpEta = ( allFz .* allFzbar ) * weight;
            interpFzBar = interpEta ./ interpFz;

            % 3) integrate fz -> Phi, fzbar -> Psi by collecting edge
            % differences.
            edgeDifferencesFz = edgeVectors .* (0.5 * (interpFz(endNodes(:,1)) + interpFz(endNodes(:,2))));
            sparseDifferencesFz = sparse([endNodes(:,1); endNodes(:,2)], ...
                                        [endNodes(:,2); endNodes(:,1)], ...
                                        [-1 * edgeDifferencesFz; edgeDifferencesFz]);
            edgeDifferencesFzBar = edgeVectors .* (0.5 * (interpFzBar(endNodes(:,1)) + interpFzBar(endNodes(:,2))));
            sparseDifferencesFzBar = sparse([endNodes(:,1); endNodes(:,2)], ...
                                        [endNodes(:,2); endNodes(:,1)], ...
                                        [-1 * edgeDifferencesFzBar; edgeDifferencesFzBar]);
            
            % Set Phi, Psi anchor values as defined in BDHI                        
            Phi = zeros( size(C,1), 1 );
            Psi = zeros( size(Phi) );
            
            mixedF = allVertices * weight;
            Phi(anchorIndex) = mixedF(anchorIndex);
            
            Psi(1) = 0;
            
            % Traverse the graph, accumulating edge values in Phi, Psi
            for i = 1:length(dfsEdges)
                currIndex = dfsEdges(i,2);
                prevIndex = dfsEdges(i,1);
                Phi(currIndex) = sparseDifferencesFz(prevIndex, currIndex) + Phi(prevIndex);
                Psi(currIndex) = sparseDifferencesFzBar(prevIndex, currIndex) + Psi(currIndex);
            end
            
            % 4) sum them together
            interpF = Phi + Psi;
        end
        
        % 5) display for various times t
        n = 100;
        increment = 1.0 / n;
        oldF = allVertices(:,1);
        for t = 0 : increment : 1
            if exist('movementVectorsHandle', 'var') > 0
                delete(movementVectorsHandle);
            end
            if t > 0 && showMovementVectors
                hold on;
                movementVectorsHandle = plotMovementVectors( ...
                    horzcat(real(oldF), imag(oldF)), ...
                    horzcat(real(newF), imag(newF)), ...
                    3 ...
                );
                pause(0);
                hold off;
                oldF = newF;
            end
            newF = interpolate(t, @quadraticSplineWeight);
            display(t);
            set(g_Deform(gid).tsh, 'Vertices', [real(newF), imag(newF)]);
            drawnow
        end
        if exist('movementVectorsHandle', 'var') > 0
            delete(movementVectorsHandle);
        end
        
    end

end