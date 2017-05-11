function coeffs = getHermiteCoefficients(eta, deta_dt)
%% Get the coefficients 
% using complex matrices: eta, deta_dt

    nVertices = size(eta, 1);
    nKeyframes = size(eta, 2);
    nKeyframePairs = nKeyframes - 1;
    
    B = [2 -2 1 1; -3 3 -2 -1; 0 0 1 0; 1 0 0 0];
    
    h = zeros( 4 , nVertices * nKeyframePairs );
    
    for i = 1:nKeyframePairs
        hBlock = [eta(:,i), eta(:,i+1), deta_dt(:,i), deta_dt(:,i+1)];
        h(1:4, nVertices*(i-1)+1 : nVertices*i) = transpose(hBlock);
    end
    coeffs = transpose(B * h);
end