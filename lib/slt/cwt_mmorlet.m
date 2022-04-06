function [C, AC] = cwt_mmorlet(x, fs, w0, freqset)
    %%
    Hfun = @(ome, w0, F) exp(-(sqrt(2)*pi*w0*(F - ome)./(5*ome)).^2);
    C = zeros(length(freqset), length(x));
    AC = zeros(length(freqset), 1);
    
    %%
    ome = freqset;
    df = (fs/length(x));
    F = (0:(length(x) - 1)) .* df;
    
    %%
    if isa(x, 'gpuArray')
        F = gpuArray(cast(F, 'single'));
        w0 = gpuArray(cast(w0, 'single'));
        ome = gpuArray(cast(ome, 'single'));
        AC = gpuArray(cast(AC, 'single'));
        C = gpuArray(cast(C, 'single'));
    end
    
    %%
    X = fft(x);
    const = w0*2*sqrt(pi)./(5.*ome);
    
    for i=1:length(ome)
        H = Hfun(ome(i), w0, F);
        
        AC(i) = sum(abs(H).^2) * const(i) * df;
        
        if abs(1 - AC(i)) < 1e-4
            psy = ifft(X .* H);
            C(i, :) = psy;
        else
            C(i, :) = NaN;
        end
    end
end