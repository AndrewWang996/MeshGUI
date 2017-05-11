function drawAnimation(meshname, trisurf_handle, numTimesPerInterval)
    info = getInformation(meshname);

    interpF = interpolateF(info, numTimesPerInterval);
    x = real(interpF);
    y = imag(interpF);
    
    xlim([ min(x(:)), max(x(:)) ]);
    ylim([ min(y(:)), max(y(:)) ]);
    
    nFrames = size(interpF, 2);
    for i = 1:nFrames
        x_vals = x(:,i);
        y_vals = y(:,i);
        set(trisurf_handle, 'Vertices', [x_vals y_vals]);
        drawnow;
    end
end

