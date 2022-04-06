function [P, H, t_P, x, freqset, dec_n, fs_x] = f0_engine(audiofilepath)
    %%
    conf = config_f0();
    addpath('./lib/cfp/');

    %%
    [s, fs_s] = audioread(audiofilepath);
    
    if size(s, 2) == 2
        s = mean(s, 2);
    end
    
    %%
    fs_x = conf.fs_ds;
    x = helper.h_downsampling(s, fs_s, fs_x);

    %%
    x = x(:)';
    t_x = (0:(numel(x) - 1))./fs_x;

    if conf.gpumode
        x = gpuArray(cast(x, 'single'));
    end

    addpath('./lib/slt/');
    [P, freqset, t_P] = h_fslt(x, t_x, fs_x, config_slt(), conf.t_dec);
    
    %%
    conf_cfp = config_cfp();
    H = cfprep(2.*P, freqset, conf_cfp.gam_X, conf_cfp.gam_Y, conf_cfp.sharpness);

    if isa(H, 'gpuArray')
        H = gather(H);
    end
    
    %%
    dec_n = round(conf.t_dec * fs_x);

    %%
    if isa(x, 'gpuArray')
        x = gather(x);
    end
end