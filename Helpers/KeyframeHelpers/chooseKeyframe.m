function whichKeyframe = ChooseKeyframe(meshname)
%     if nargin < 1
%         meshname = 'square';
%         disp(strcat('please provide a meshname to ChooseKeyframe function.',...
%             ' Defaulting to square'));
%     end

    f = figure;
    showKeyframeSlider = uicontrol(f,'Style','slider',...
        'String','Drag to choose keyframe.',...
        'Position',[0 0 120 20],...
        'Callback', @DiscretizeSliderMovement);
    addlistener(showKeyframeSlider, ...
        'ContinuousValueChange', @KeyframePicker);
    
    doneButton = uicontrol(f,...
        'String','Finish',...
        'Position',[200 0 120 20],...
        'Callback', @FinishKeyframeSelection);

    
    keyframe_1 = getKeyframe(meshname, 1);
    V_1 = keyframe_1.Vertices;
    F_1 = keyframe_1.Faces;
    h = trisurf(F_1,V_1(:,1),V_1(:,2),zeros(size(V_1,1),1), ...
        'FaceColor','interp');
    view(2);
    axis equal;
    axis manual;
    
    numKeyframes = countKeyframes(meshname);
    
    function FinishKeyframeSelection(hObject, eventdata, handles)
        if numKeyframes == 0
            disp('no keyframes found');
            return
        end
        val = capValue( showKeyframeSlider.Value );
        whichKeyframe = floor(val * numKeyframes) + 1;
        close(gcf);
    end

    function KeyframePicker(src, event)
        if numKeyframes == 0
            disp('no keyframes found');
            return 
        end
        position = src.Value;
        i = round( (numKeyframes - 1) * position + 1);
        keyframe = getKeyframe(meshname, i);
        set(h, 'Vertices', keyframe.Vertices);
    end

    function DiscretizeSliderMovement(hObject, eventdata, handles)
        if numKeyframes == 0
            disp('no keyframes found');
            return
        end
        val = capValue(hObject.Value);
        
        i = floor(val * numKeyframes) + 1;
        roundedVal = (i - 1) / (numKeyframes - 1);
        hObject.Value = roundedVal;
        disp(strcat('keyframe:' , int2str(i)));
    end

    uiwait(f);
end

    
function value = capValue(val)
    eps = 1e-6;
    value = val;
    if val < 0
        value = 0;
    elseif val > 1 - eps
        value = 1 - eps;
    end
end