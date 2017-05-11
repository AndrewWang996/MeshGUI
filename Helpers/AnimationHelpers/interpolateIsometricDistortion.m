function isom_dist = interpolateIsometricDistortion(info, numTimesPerInterval)
    logFz = info.allLogFz;
    allEta = info.allEta;
    

    interpFz = exp( getInterpolatedPointsBezier( logFz, numTimesPerInterval ) );
    interpEta = getInterpolatedPointsBezier(allEta, numTimesPerInterval);
    interpFzBar = interpEta ./ conj(interpFz);
    
    sig_a = abs(interpFz) + abs(interpFzBar);
    sig_b = abs( abs(interpFz) - abs(interpFzBar) );
    
    
    isom_dist = sig_a + (1 ./ sig_b);
end