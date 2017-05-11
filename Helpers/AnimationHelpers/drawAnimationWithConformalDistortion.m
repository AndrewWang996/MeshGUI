function drawAnimationWithConformalDistortion(meshname, trisurf_handle, numTimesPerInterval)
    info = getInformation(meshname);

    interpF = interpolateF(info, numTimesPerInterval);
    conformal_distortion = interpolateConformalDistortion(info, numTimesPerInterval);
    
    x = real(interpF);
    y = imag(interpF);
    
    xlim([ min(x(:)), max(x(:)) ]);
    ylim([ min(y(:)), max(y(:)) ]);
    
    nFrames = size(interpF, 2);
    for i = 1:nFrames
        x_vals = x(:,i);
        y_vals = y(:,i);
        conf_dist = conformal_distortion(:,i);
        set(trisurf_handle, 'Vertices', [x_vals y_vals]);
        set(trisurf_handle, 'FaceVertexCData', conf_dist);
        drawnow;
    end
end

