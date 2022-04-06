function C = fslt(x, fs, w0, freqset, o, o_int, al)
    %%
    lnC = zeros(length(freqset), length(x));
    
    if isa(x, 'gpuArray')
        lnC = gpuArray(cast(lnC, 'single'));
    end
    
    %%
    o_max = max(o_int);
    
    f = waitbar(0, 'Please wait...');
    
    for i=1:o_max
        waitbar(i/o_max, f, 'Calculating fractional superlet of the current segment');
        
        idx = i <= o_int;
        freqset_i = freqset(idx);
        
        C_i = cwt_mmorlet(x, fs, w0*i, freqset_i);
        
        lnC(idx, :) = lnC(idx, :) + log(C_i);
        
        idx = o_int == i;
        
        if ~all(idx == 0)
            freqset_j = freqset(idx);
            C_j = cwt_mmorlet(x, fs, w0*(i + 1), freqset_j);

            al_j = al(idx)';
            lnC(idx, :) = lnC(idx, :) + bsxfun(@times, al_j, log(C_j));
        end
    end
    
    waitbar(1, f, 'Finishing');
    close(f)

    %%
    C = exp(bsxfun(@times, 1./o', lnC));
end