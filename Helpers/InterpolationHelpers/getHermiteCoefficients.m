function coeffs = getHermiteCoefficients(eta, deta_dt)
%% Get the coefficients 
% using complex matrices: eta, deta_dt

    nVertices = size(eta, 1);
    nKeyframes = size(eta, 2);
    nKeyframePairs = nKeyframes - 1;
    
    B = [2 -2 1 1; -3 3 -2 -1; 0 0 1 0; 1 0 0 0];
    repB = repmat(B, 1, nVertices);
    
    repH = zeros( 4 * nKeyframePairs, nVertices );
    repH(1:4:end, :) = eta(1:nKeyframePairs, :);
    repH(2:4:end, :) = eta(2:nKeyframePairs+1, :);
    repH(3:4:end, :) = deta_dt(1:nKeyframePairs, :);
    repH(4:4:end, :) = deta_dt(2:nKeyframePairs+1, :);
    
    coeffs = repB * repH;
end