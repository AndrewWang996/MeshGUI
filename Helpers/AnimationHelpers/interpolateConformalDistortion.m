function conf_dist = interpolateConformalDistortion(info, numTimesPerInterval)
    logFz = info.allLogFz;
    allEta = info.allEta;
    

    interpFz = exp( getInterpolatedPointsBezier( logFz, numTimesPerInterval ) );
    interpEta = getInterpolatedPointsBezier(allEta, numTimesPerInterval);
    interpFzBar = interpEta ./ conj(interpFz);
    
    sig_a = abs(interpFz) + abs(interpFzBar);
    sig_b = abs( abs(interpFz) - abs(interpFzBar) );
    
    K = sig_a ./ sig_b;
    k = (K - 1) ./ (K + 1);
    
    conf_dist = k;
end