function demo_fslt2
    rng(15);

    %%
    fs = 5000;
    T = 10;
    K = randi(10, 1);
    
    [x, t, f0, ~, A] = testlocsignal(fs, T, K);
    
    %%
    w0 = 4;
    gpumode = true;
    freqset = 440 .* 2.^((-530:300)./120);
    
    o_min = 16;
    o_max = 50;

    %%
    o = o_min + (o_max - o_min) * (freqset - freqset(1))./(freqset(end) - freqset(1));
    o_int = floor(o);
    al = o - o_int;
    
    lnC = zeros(length(freqset), length(x));
    if gpumode
        lnC = gpuArray(cast(lnC, 'single'));
        al = gpuArray(cast(al, 'single'));
        o = gpuArray(cast(o, 'single'));
    end
    
    for i=1:o_max
        idx = i <= o_int;
        freqset_i = freqset(idx);
        
        C_i = cwt_mmorlet(x, fs, w0*i, freqset_i, gpumode);
        
        lnC(idx, :) = lnC(idx, :) + log(C_i);
        
        idx = o_int == i;
        
        if ~all(idx == 0)
            freqset_j = freqset(idx);
            C_i = cwt_mmorlet(x, fs, w0*(i + 1), freqset_j, gpumode);

            al_j = al(idx)';
            lnC(idx, :) = lnC(idx, :) + bsxfun(@times, al_j, log(C_i));
        end
    end
    
    C = exp(bsxfun(@times, 1./o', lnC));

    P = abs(C).^2;
    
    %%
    dec_t = 0.01;
    dec_n = round(dec_t * fs);
    
    %%
    freqset_B = freqset;
    B = P .* 0;
    
    for k=1:K
        [~, idx] = min(abs(f0(k) - freqset));
        freqset_B(idx) = f0(k);
        B(idx, :) = A(k, :);
    end
    
    %%
    figure(3);
    subplot(2, 1, 1);
    surf(t(1:dec_n:end), freqset, 2.*sqrt(P(:, 1:dec_n:end)), 'edgecolor', 'none');
    view(0, 90);
    axis tight;
    ylim([f0(1)*0.95 f0(end)*1.05]);
    set(gca, 'YScale', 'log');
    
    subplot(2, 1, 2);
    surf(t(1:dec_n:end), freqset_B, B(:, 1:dec_n:end), 'edgecolor', 'none');
    view(0, 90);
    axis tight;
    ylim([f0(1)*0.95 f0(end)*1.05]);
    set(gca, 'YScale', 'log');
    title(sprintf('K = %d', K));
end