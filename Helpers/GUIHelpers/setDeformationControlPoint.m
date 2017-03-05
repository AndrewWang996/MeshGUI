function setDeformationControlPoint(src, event)
    global deforming
    global deformStartDefined
    if exist('deforming', 'var') ~= 1
        return
    elseif deforming ~= 1
        return
    end
    display('set deformation control point')
    clickpoint = get(gca,'currentpoint');

    x = clickpoint(1,1,1);
    y = clickpoint(1,2,1);

    hold on;
    plot(x,y,'*');
    hold off;

    deforming = true;
    deformStartDefined = true;

end