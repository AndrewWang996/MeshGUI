
addpath Helpers;
addpath Helpers/KeyframeHelpers;
addpath Helpers/PlotHelpers;
addpath Helpers/GUIHelpers;
addpath Helpers/OptimizationHelpers;
addpath Helpers/GraphHelpers;
addpath Helpers/InterpolationHelpers;
addpath Helpers/InterpolationHelpers/BezierHelpers;
addpath Helpers/InterpolationHelpers/HermiteHelpers;
addpath Helpers/AnimationHelpers;


addpath WeightFunctions;

addpath Meshes;

axis equal;



meshname = 'red_dragon';
[V,F] = getMesh(meshname);
cagepts = getCage(meshname);

set(gcf, 'UserData', struct('meshname', meshname));
set(gcf, 'Position', [0 200 900 500]);

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

showConformalDistortion = uicontrol(gcf,'Style','pushbutton',...
    'String','Conf. Dist.',...
    'Position', [500 0 90 20], ...
    'Callback', @ShowConformalDistortion);

showIsometricDistortion = uicontrol(gcf,'Style','pushbutton',...
    'String','Isom. Dist.',...
    'Position', [600 0 90 20], ...
    'Callback', @ShowIsometricDistortion);

showKeyframeButton = uicontrol(gcf,'Style','slider',...
    'String','Save Keyframe',...
    'Position',[700 0 120 20]);
addlistener(showKeyframeButton, ...
    'ContinuousValueChange', @ShowKeyframe);

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
        anchorPositions = complex(vertices(anchorIndices, 1), vertices(anchorIndices, 2));

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
    
    function ShowConformalDistortion(src,event)
        drawAnimationWithConformalDistortion(meshname, g_Deform(gid).tsh, 100);
    end

    function ShowIsometricDistortion(src,event)
        drawAnimationWithIsometricDistortion(meshname, g_Deform(gid).tsh, 100);
    end
    
    function SetVelocity(src,event)
        disp( strcat('Select pairs of points. First point a vertex on mesh.', ...
            ' Second point is endpoint of velocity vector.') );
        
        whichKeyframe = chooseKeyframe(meshname);
        keyframe = getKeyframe(meshname, whichKeyframe);
        vertices = keyframe.Vertices(:,1:2);
        faces = keyframe.Faces;
        
        figure
        trisurf(faces, vertices(:,1), vertices(:,2), ...
            zeros(size(vertices,1),1));
        view(2);
        axis equal;
        axis manual;
        
        while 1
            [ptsX, ptsY] = getpts; 
            if mod(size(ptsY, 1), 2) ~= 0
                disp( strcat('please remember to select an ', ...
                    'even number of points.\n') );
                continue;
            end
            break;
        end

        complexPts = complex(ptsX, ptsY);

        ptsFrom = complexPts(1:2:length(complexPts), :);
        ptsTo = complexPts(2:2:length(complexPts), :);

        indices = getIndex(ptsFrom, vertices);
        hold on;
        df_dt_handles = plotMovementVectors( ...
            vertices(indices,1:2), ...
            horzcat(real(ptsTo), imag(ptsTo)), ...
            1 ...
        );
        
        
        velocities_unamplified = ptsTo - complex(vertices(indices,1), vertices(indices,2));
        velocities = 50 * velocities_unamplified;
        % use the first vector as anchor
        anchorIndex = indices(1);
        indices = indices(2:end);
        velocities = velocities - velocities(1,:);
        velocities = velocities(2:end,:);
        [H,Aeq,beq] = poseOptimizationProblem(meshname, indices, velocities, whichKeyframe, anchorIndex);
        deta_dt = solveOptimizationProblem(H,Aeq,beq);
        
        save_deta_dt(meshname, whichKeyframe, deta_dt);
        
        
        deta_dt_handles = plotMovementVectors( ...
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
        drawAnimation(meshname, g_Deform(gid).tsh, 100);
    end

end