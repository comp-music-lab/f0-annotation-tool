function [P, freqset, t_P] = h_fslt(x, t_x, fs_x, conf, dec_t)
    %%
    w0 = conf.w0;
    o_min = conf.o_min;
    o_max = conf.o_max;
    freqset = conf.freqset;
    N_max = conf.N_max;
    
    %%
    if strcmp(conf.order, 'i')
        [o, o_int, al] = iadaptiveorder(freqset, o_min, o_max);
    else
        [o, o_int, al] = adaptiveorder(freqset, o_min, o_max);
    end
    
    %%
    M = floor(numel(x)*numel(freqset)/N_max);
    k = round(linspace(1, numel(freqset), 2 + M));

    %%
    dec_n = round(fs_x * dec_t);
    P = zeros(numel(freqset), numel(x(1:dec_n:end)));
    
    idx_st = k(1);

    for i=2:numel(k)
        idx_ed = k(i);
        I = idx_st:idx_ed;
        C = fslt(x, fs_x, w0, freqset(I), o(I), o_int(I), al(I));

        P_i = gather(2.*abs(C).^2);
        P(I, :) = P_i(:, 1:dec_n:end);

        idx_st = idx_ed + 1;
    end

    %%
    t_P = t_x(1:dec_n:end);
end