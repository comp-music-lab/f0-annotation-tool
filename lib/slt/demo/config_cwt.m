function cwtconf = config_cwt
    freqset = 440 .* 2.^((-530:500)./120);
    
    name = 'wavelet';
    w0 = 24;
    o_min = 1;
    o_max = 1;
    
    cwtconf = struct('name', name, 'freqset', freqset, 'w0', w0, 'o_min', o_min, 'o_max', o_max);
end