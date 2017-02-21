function setDeformationControlEnd(src, event)
    global deforming
    global deformStartDefined
    if exist('deforming', 'var') ~= 1 || ~ deforming
        return
    end

    display('set deformation control end')

    clickpoint = get(gca,'currentpoint');
    x = clickpoint(1,1,1);
    y = clickpoint(1,2,1);

    hold on;
    plot(x,y,'*');
    hold off;



    deforming = false;
    deformStartDefined = false;

end