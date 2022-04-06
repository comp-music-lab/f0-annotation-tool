function demo_modifiedMorlet3
    rng(15);

    %%
    fs = 5000;
    T = 10;
    K = randi(10, 1);
    
    [x, t, f0, ~, A] = testlocsignal(fs, T, K);
    
    %%
    w0 = 32;
    gpumode = false;
    freqset = 440 .* 2.^((-530:300)./120);
    
    [C, AC] = cwt_mmorlet(x, fs, w0, freqset, gpumode);
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
    figure(1);
    plot(AC);
    
    figure(2);
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