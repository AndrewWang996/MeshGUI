function isom_dist = interpolateIsometricDistortion(info, numTimesPerInterval)
    logFz = info.allLogFz;
    allEta = info.allEta;
    all_deta_dt = info.all_deta_dt;
    

    interpFz = exp( getInterpolatedPointsBezier( logFz, numTimesPerInterval ) );
    interpEta = getInterpolatedPointsHermite(allEta, all_deta_dt, numTimesPerInterval);
    interpFzBar = interpEta ./ conj(interpFz);
    
    sig_a = abs(interpFz) + abs(interpFzBar);
    sig_b = abs( abs(interpFz) - abs(interpFzBar) );
    
    
    isom_dist = sig_a + (1 ./ sig_b);
end