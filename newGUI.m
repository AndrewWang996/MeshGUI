
addpath Helpers;
addpath Helpers/KeyframeHelpers;
addpath Helpers/PlotHelpers;
addpath Helpers/GUIHelpers;

addpath WeightFunctions;

addpath Meshes;

axis equal;

%{
fprintf( ...
    ['\nCLICK on mesh at each location where you would like to add a ' ...
    'point handle.\n' ...
    'Press ENTER when finished.\n\n']);
%}


meshname = 'simple';
[V,F] = getMesh(meshname);
cagepts = getCage(meshname);

set(gcf, 'UserData', struct('meshname', meshname));

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
    % plot bones
    hold on;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % cache fz, fzbar for the identity function f(z) = z
    g_Deform(gid).fz = ones( length(V), 1 );
    g_Deform(gid).fzbar = zeros( length(V), 1 );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveKeyframeButton = uicontrol(gcf,'Style','pushbutton',...
    'String','Save Keyframe',...
    'Position',[50 0 90 20],...
    'Callback', @SaveKeyframe);

showKeyframeButton = uicontrol(gcf,'Style','slider',...
    'String','Save Keyframe',...
    'Position',[400 0 120 20]);
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

    function StartDeformation(src,event)
        display('click pairs of points, 1st on the shape, 2nd on the desired new location')

        [vertices, faces] = getMesh(meshname);
        cagePts = getCage(meshname);

        figure
        trimesh(faces, vertices(:,1), vertices(:,2));
        [ptsX, ptsY] = getpts;
        close(gcf)

        complexPts = ptsX + 1i * ptsY;

        ptsFrom = complexPts(1:2:length(complexPts), :);
        ptsTo = complexPts(2:2:length(complexPts), :);

        indices = getIndex(ptsFrom, vertices);
        anchorIndices = getAnchorIndices(meshname);
        anchorPositions = vertices(anchorIndices, 1) + 1i * vertices(anchorIndices, 2);

        [newVerticesComplex, fz, fzbar] = deformBoundedDistortion([indices; anchorIndices], [ptsTo; anchorPositions], vertices, faces, cagePts);
        newVertices = [real(newVerticesComplex), imag(newVerticesComplex)];

        set(g_Deform(gid).tsh, 'Vertices', newVertices);
        g_Deform(gid).fz = fz;
        g_Deform(gid).fzbar = fzbar;
    end

    hold off;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up interaction variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % keep track of window xmin, xmax, ymin, ymax
    win_min = min([C(:,1:2); V(:,1:2)]);
    win_max = max([C(:,1:2); V(:,1:2)]);
    % keep track of down position
    down_pos = [];
    % keep track of last two drag positions
    drag_pos = [];
    last_drag_pos = [];
    % keep track of mesh vertices at mouse down
    down_V = [];
    % keep track of index of selected control point
    ci = [];
    % type of click ('left','right')
    down_type  = '';

    fprintf( ...
        ['DRAG a control point to deform the mesh.\n' ...
        'RIGHT CLICK DRAG a control point to rotate point handles.\n\n']);

    return

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Callback functions for keyboard and mouse
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Callback for mouse down on control points
    function oncontrolsdown(src,event)
        % get current mouse position, and remember old one
        down_pos=get(gca,'currentpoint');
        down_pos=[down_pos(1,1,1) down_pos(1,2,1)];
        last_drag_pos=down_pos;
        drag_pos=down_pos;
        % keep track of control point positions at mouse down
        g_Deform(gid).new_C = [get(C_plot,'XData')' get(C_plot,'YData')'];
        % get index of closest control point
        [~,ci] =  ...
            min(sum((g_Deform(gid).new_C(:,1:2) - ...
            repmat(down_pos,size(g_Deform(gid).new_C,1),1)).^2,2));
        % keep track of mesh vertices at mouse down
        down_V = get(g_Deform(gid).tsh,'Vertices');
        down_V = down_V(:,1:2);
        
        % tell window that drag and up events should be handled by controls
        set(gcf,'windowbuttonmotionfcn',@oncontrolsdrag)
        set(gcf,'windowbuttonupfcn',@oncontrolsup)
        set(gcf,'KeyPressFcn',@onkeypress)
        if(strcmp('normal',get(gcf,'SelectionType')))
            % left-click
            down_type = 'left';
        else
            % other (right) click
            down_type = 'right';
        end
        
    end

% Callback for mouse drag on control points
    function oncontrolsdrag(src,event)
        % keep last drag position
        last_drag_pos=drag_pos;
        % get current mouse position
        drag_pos=get(gca,'currentpoint');
        drag_pos=[drag_pos(1,1,1) drag_pos(1,2,1)];
        if(strcmp('left',down_type))
            % move selected control point by drag offset
            g_Deform(gid).new_C(ci,:) = ...
                g_Deform(gid).new_C(ci,:) + drag_pos-last_drag_pos;
        end
        update_positions();
    end


    function update_positions()
        % update display positions
        set(C_plot,'XData',g_Deform(gid).new_C(:,1));
        set(C_plot,'YData',g_Deform(gid).new_C(:,2));
        
        % update mesh positions
        new_V = g_Deform(gid).tsh.Vertices;
        for i = 1:length(g_Deform(gid).indices)
            index = g_Deform(gid).indices(i);
            new_V(index,1) = g_Deform(gid).new_C(i,1);
            new_V(index,2) = g_Deform(gid).new_C(i,2);
        end
        
        % update mesh positions
        set(g_Deform(gid).tsh,'Vertices',new_V(:,1:2));
    end

% Callback for mouse release of control points
    function oncontrolsup(src,event)
        % Tell window to handle drag and up events itself
        set(gcf,'windowbuttonmotionfcn','');
        set(gcf,'windowbuttonupfcn','');
        cur_V = get(g_Deform(gid).tsh,'Vertices');
        cur_V = cur_V(:,1:2);
        
        % scale window to fit
        win_min = min([win_min; cur_V]);
        win_max = max([win_max; cur_V]);
        axis(reshape([win_min;win_max],1,2*size(cur_V,2)))
    end

    function onkeypress(src,event)
        if(strcmp(event.Character,'r'))
            
            % Refresh the state of the mesh
            g_Deform(gid).new_C = C;
            update_positions();
            
        elseif(strcmp(event.Character,'u'))
            update_positions();
        end
    end


    function SaveKeyframe(src, event)
        copyOfCurrent = {};
        copyOfCurrent.Vertices = g_Deform(gid).tsh.Vertices;
        copyOfCurrent.Faces = g_Deform(gid).tsh.Faces;
        copyOfCurrent.fz = g_Deform(gid).fz;
        copyOfCurrent.fzbar = g_Deform(gid).fzbar;
        
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
%         set(C_plot,'XData',keyframe.Vertices(:,1));
%         set(C_plot,'YData',keyframe.Vertices(:,2));
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
        [endNodes, weights, predecessor] = getSpanningTree(meshname);
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
                movementVectorsHandle = plotMovementVectors(oldF, newF, 3);
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