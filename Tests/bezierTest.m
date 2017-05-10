function bezierTest()
    addpath Helpers/InterpolationHelpers/BezierHelpers;
    addpath Helpers/InterpolationHelpers/HermiteHelpers;
    
    v1 = [0 + 0i; 0 + 1i; 1 + 0i];
    v2 = [0 + 2i; 1 + 2i; 0 + 1i];
    v3 = [2 + 2i; 2 + 1i; 1 + 2i];
    v4 = [2 + 0i; 1 + 0i; 2 + 1i];
    
    interpolated_points = horzcat(v1, v2, v3, v4);
    
    numTimesPerInterval = 1;
    
    points = getInterpolatedPointsBezier(interpolated_points, numTimesPerInterval);
    
    figure;
    hold on;
    disp([ min( min( real(points) ) ) , max( max( real(points) ) ) ]);
    xlim([ min( min( real(points) ) ) , max( max( real(points) ) ) ]);
    ylim([ min( min( imag(points) ) ) , max( max( imag(points) ) ) ]);
    for t = 1:size(points,2)
        points_snapshot = points(:,t);
        x = real(points_snapshot);
        y = imag(points_snapshot);
        plot( x , y );
        drawnow;
    end
    hold off;
    
end