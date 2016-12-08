function pts = boundaryPicker(filepath)
    if nargin < 1
        [V,F] = readOff('Meshes/elephant.off');
    else
        [V,F] = readOff(filepath);
    end
	trimesh(F, V(:,1), V(:,2), 'color', 'k')
    
    fprintf('select a boundary for the mesh in counterclockwise order and then press enter');
    
    [ptsX, ptsY] = getpts; % gets the points in the order that they are selected
    pts = ptsX + 1i * ptsY;
    
    % plot(ptsX, ptsY);
    graphComplex(pts)

end


