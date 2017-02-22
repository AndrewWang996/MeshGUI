function indices = anchorPicker(filepath)
    if nargin < 1
        [V,F] = readOff('Meshes/elephant/elephant.off');
    else
        [V,F] = readOff(filepath);
    end
	trimesh(F, V(:,1), V(:,2), 'color', 'k');
    
    fprintf('select 2 anchor points in the interior\n');
    while 1
        [ptsX, ptsY] = getpts; % gets the points in the order that they are selected
        if size(ptsY, 1) ~= 2
            fprintf('please remember to select TWO anchor points...\n');
            continue;
        end
        break;
    end
    
    indices = knnsearch(V, [ptsX, ptsY, zeros(size(ptsY, 1), 1)]);

    % plot(ptsX, ptsY);
    hold on;
    scatter( V(indices, 1), V(indices, 2) )
    hold off;

    pause(1.000);

    close;
    
end


