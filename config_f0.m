function conf = config_f0
    if gpuDeviceCount > 0
        gpumode = true;
    else
        gpumode = false;
    end

    conf = struct(...
        'fs_ds', 16000,...
        't_dec', 0.005,...
        'gpumode', gpumode,...
        'outputdir', './output/'...
        );
end