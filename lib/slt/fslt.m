function C = fslt(x, fs, w0, freqset, var1, var2, var3)
    if nargin == 6
        o_min = var1;
        o_max = var2;

        [o, o_int, al] = adaptiveorder(freqset, o_min, o_max);
    elseif nargin == 7
        o = var1;
        o_int = var2;
        al = var3;
    end
    
    al = al(:);
    o_int_max = max(o_int);

    %%
    lnC = zeros(length(freqset), length(x));
    
    if isa(x, 'gpuArray')
        lnC = gpuArray(cast(lnC, 'single'));
    end
    
    %%
    figw = waitbar(0, 'Please wait...');
    
    for i=1:o_int_max
        waitbar(i/o_int_max, figw, 'Calculating the fractional superlet of the current segment');
        
        % select target frequency band
        idx = find(i <= o_int);
        freqset_i = freqset(idx);
        
        C_i = cwt_mmorlet(x, fs, w0*i, freqset_i);
        lnC(idx, :) = lnC(idx, :) + log(C_i);
        
        % geometric mean
        idx = o_int == i;
        
        if ~all(idx == 0)
            freqset_j = freqset(idx);
            al_j = al(idx);

            if ~all(al_j == 0)
                C_j = cwt_mmorlet(x, fs, w0*(i + 1), freqset_j);
                lnC(idx, :) = lnC(idx, :) + bsxfun(@times, al_j, log(C_j));
            end
        end
    end
    
    %%
    waitbar(1, figw, 'Finishing');
    close(figw)
    
    C = exp(bsxfun(@times, 1./o, lnC));
end