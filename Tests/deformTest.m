function testDeform()
    [mesh,faces] = getMesh('red_dragon');
    cage = getCage('red_dragon');
    graphComplex(cage);
    
    indices = [100; 37; 300; 600; 999; 1200; 1450];
    points = mesh(indices , :);
    
    hold on;
    for i = 1:size(points,1)
        x = points(i,1);
        y = points(i,2);
        plot(x,y, '');
        text(x,y,num2str(i));
    end
    hold off;
    
    complexPoints = points(:,1) + 1i * points(:,2);
    complexPoints(1) = 0 + 1i * 1;
    
    [newVertices, fz, fzbar, phi, psi] = deformBoundedDistortion(indices, complexPoints, mesh, faces, cage);
    
    figure
    trimesh(faces,real(newVertices),imag(newVertices),'color','k');
    drawnow;
end
