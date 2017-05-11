function getControlPointsTest()
    interpolated_points = [0,0; 1,2; 2,10; 5,6];
    control_points = getControlPoints(interpolated_points);
    
    hold on;
    plot(interpolated_points(:,1), interpolated_points(:,2));
    plot(control_points(:,1), control_points(:,2));
    hold off;
end