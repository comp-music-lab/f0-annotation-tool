function demo_modifiedMorlet2
    %%
    fs = 1000;
    dt = 1/fs;
    T = 5;
    t = 0:dt:T;
    K = randi(10, 1);
    
    A = rand(K, 1);
    f = rand(K, 1) .* (fs/2 * 0.75);
    f = sort(f);
    
    x = zeros(1, length(t));
    for k=1:K
        x = x + A(k).*sin(2*pi*f(k).*t + 2*pi*rand);
    end
    
    %%
    w0 = 32;
    gpumode = false;
    freqset = f;
    
    [C, AC] = cwt_mmorlet(x, fs, w0, freqset, gpumode);
    P = abs(C).^2;
    
    %%
    figure(1);
    plot(AC);
    
    figure(1);
    for k=1:K
        subplot(K, 1, k);
        plot(t, 2.*sqrt(P(k, :))); hold on;
        plot([t(1) t(end)], [A(k) A(k)], '-.m'); hold off;
        title(sprintf('f = %3.3f', f(k)));
    end
end