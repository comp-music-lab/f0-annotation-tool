function cfpconf = config_cfplog
    gam_X = 0.01;
    gam_Y = 1.1;
    sharpness = 32;
    
    cfpconf = struct('gam_X', gam_X, 'gam_Y', gam_Y, 'sharpness', sharpness);
end