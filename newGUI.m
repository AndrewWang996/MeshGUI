
[x,y] = meshgrid(1:10, 1:10);
tri = delaunay(x,y);
handle = trimesh(tri,x,y);

axis equal;
fprintf( ...
    ['\nCLICK on mesh at each location where you would like to add a ' ...
    'point handle.\n' ...
    'Press ENTER when finished.\n\n']);


global keyframes;
keyframes = {};

V = [reshape(x, [numel(x),1]), reshape(y, [numel(y),1])];
F = tri;
simple_deform(V, F)




function simple_deform(varargin)

V = varargin{1};
% face indices of mesh
F = varargin{2};
% control vertices of skeleton
C = V;
% what indices of V do C correspond to
I = 1:length(V);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set default point handles
P = 1:size(C,1);
% Be sure that control vertices are in 2D
if(size(C,2) == 3)
    C = C(:,1:2);
end


% number of point handles
np = numel(P);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global g_Deform;
gid = numel(g_Deform)+1;
% store indices I
g_Deform(gid).indices = I;
% keep track of control positions at mouse down
g_Deform(gid).new_C = [];
% keep track of rotations stored at each control point, for 2D this is a m
% by 1 list of angles
g_Deform(gid).R = zeros(np,1);
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

% plot the control points (use 3D plot and fake a depth offset by pushing
% control points up in z-direction)

C_plot = scatter3( ...
    C(:,1),C(:,2),0.1+0*C(:,1), ...
    'o','MarkerFaceColor',[0.9 0.8 0.1], 'MarkerEdgeColor','k',...
    'LineWidth',2,'SizeData',20, ...
    'ButtonDownFcn',@oncontrolsdown);
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
        
        % tell window that drag and up eventents should be handled by controls
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
        else
            [found, iP] = ismember(ci,P);
            if(found)
                g_Deform(gid).R(iP) = ...
                    g_Deform(gid).R(iP) + 2*pi*(drag_pos(1)-last_drag_pos(1))/100;
            end
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
        % Tell window to handle drag and up eventents itself
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
            g_Deform(gid).R = zeros(np,1);
            update_positions();
            
        elseif(strcmp(event.Character,'u'))
            update_positions();
        end
    end


    function SaveKeyframe(src, event)
        global keyframes;
        copyOfCurrent = {};
        copyOfCurrent.Vertices = g_Deform(gid).tsh.Vertices;
        copyOfCurrent.Faces = g_Deform(gid).tsh.Faces;
        
        keyframes = [keyframes, copyOfCurrent];
        fprintf('Saved 1 new keyframe. %d total keyframes.\n', numel(keyframes));
    end

    function ShowKeyframe(src,event)
        global keyframes;
        numKeyframes = numel(keyframes);
        if numKeyframes == 0
            return 
        end
        position = src.Value;
        whichKeyframe = round( (numKeyframes - 1) * position + 1);
        % display(whichKeyframe);
        set(g_Deform(gid).tsh,'Vertices',keyframes{whichKeyframe}.Vertices(:,1:2));
        set(C_plot,'XData',keyframes{whichKeyframe}.Vertices(:,1));
        set(C_plot,'YData',keyframes{whichKeyframe}.Vertices(:,2));
    end
    
    function ShowAnimation(src,event)
        global keyframes;
        numKeyframes = numel(keyframes);
        allVertices = zeros(size(C,1), numKeyframes);
       
        for whichKeyframe = 1:numKeyframes
            allVertices(:,whichKeyframe) = complex(...
                keyframes{whichKeyframe}.Vertices(:,1),...
                keyframes{whichKeyframe}.Vertices(:,2)...
            );
        end
        
        display(allVertices);
        n = 10;
        
        for t = 0:1.0/n:1
            w = linearWeight(2,t);
            new_Vertices = allVertices * w;
            xdata = real(new_Vertices);
            ydata = imag(new_Vertices);
            set(C_plot,'XData',xdata);
            set(C_plot,'YData',ydata);
            set(g_Deform(gid).tsh, 'Vertices', [xdata,ydata]);
            pause(0.1);
        end
        
    end

end