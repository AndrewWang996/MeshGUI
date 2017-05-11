function conf_dist = interpolateConformalDistortion(info, numTimesPerInterval)
    logFz = info.allLogFz;
    allEta = info.allEta;
    all_deta_dt = info.all_deta_dt;
    

    interpFz = exp( getInterpolatedPointsBezier( logFz, numTimesPerInterval ) );
    interpEta = getInterpolatedPointsHermite(allEta, all_deta_dt, numTimesPerInterval);
    interpFzBar = interpEta ./ conj(interpFz);
    
    sig_a = abs(interpFz) + abs(interpFzBar);
    sig_b = abs( abs(interpFz) - abs(interpFzBar) );
    
    K = sig_a ./ sig_b;
    k = (K - 1) ./ (K + 1);
    
    conf_dist = k;
end