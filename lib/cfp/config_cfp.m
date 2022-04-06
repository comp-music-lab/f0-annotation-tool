function cfpconf = config_cfp
    gam_X = 0.1;
    gam_Y = 1.1;
    sharpness = 8;
    
    cfpconf = struct('gam_X', gam_X, 'gam_Y', gam_Y, 'sharpness', sharpness);
end