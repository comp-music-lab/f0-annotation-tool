function [C, freqset] = fslt_wrapper(x, fs, wconf, maxwid, gpumode)
    %%
    if gpumode
        x = gpuArray(cast(x, 'single'));
    end
    
    %%
    T = floor(6.5 * (22050/fs));
    N = round(T*fs);
    
    %%
    n_start = 1;
    n_end = n_start + N - 1;
    overlap = round(fs*(maxwid + 0.1));
    L = length(x);
    
    %%
    w0 = wconf.w0;
    o_max = wconf.o_max;
    o_min = wconf.o_min;
    freqset = wconf.freqset;
    freqset(freqset > (fs/2)) = [];
    
    [o, o_int, al] = adaptiveorder(freqset, o_min, o_max);
    
    %%
    C = [];
    
    f = waitbar(0, 'Please wait...');
    
    while true
        waitbar(n_end/L, f, 'Running fractional superlet per segment');
        
        x_i = x(n_start:n_end);
        
        C_i = fslt(x_i, fs, w0, freqset, o, o_int, al);
        
        if (n_end - n_start + 1) == N
            C = [C gather(C_i(:, 1:end - overlap))];
        else
            C = [C gather(C_i)];
            break;
        end
        
        n_start = n_end - overlap + 1;
        n_end = n_start + N - 1;
        
        if n_end > L
            n_end = L;
        end
    end
    
    waitbar(1, f, 'Finishing');
    close(f)
end