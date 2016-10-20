[x,y] = meshgrid(1:2, 1:2);
tri = delaunay(x,y);
handle = trimesh(tri,x,y);

line = handle(1);
display(line)

line.XData = [2, 0, 2, 2];
display(line)

axis equal;
fprintf( ...
  ['\nCLICK on mesh at each location where you would like to add a ' ...
  'point handle.\n' ...
  'Press ENTER when finished.\n\n']);
% User clicks many times on mesh at locations of control points
try
  [Cx,Cy] = getpts;
catch e
  % quit early, stop script
  return
end


line.LineStyle = '--';
line.ButtonDownFcn = @ButtonDownFunction;

c = uicontextmenu;
line.UIContextMenu = c;
get(line)

m1 = uimenu(c,'Label','dashed','Callback',@setlinestyle);
m2 = uimenu(c,'Label','dotted','Callback',@setlinestyle);
m3 = uimenu(c,'Label','solid','Callback',@setlinestyle);

function setlinestyle(source,eventdata)
    disp(eventdata)
    line = eventdata.Source.Parent.Parent.CurrentObject;
%     disp(line)
    switch source.Label
        case 'dashed'
            line.LineStyle = '--';
        case 'dotted'
            line.LineStyle = ':';
        case 'solid'
            line.LineStyle = '-';
    end
end

function ButtonDownFunction(hObject, eventdata)
% hObject handle to axes1 (see GCBO)
% eventdata reserved - to be defined in a future version of MATLAB

display(hObject)
display(eventdata)

end
