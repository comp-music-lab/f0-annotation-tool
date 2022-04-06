function demo_fslt1
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
    w0 = 4;
    gpumode = false;
    freqset = f;
    
    o_min = 20;
    o_max = 50;
    
    if K == 1
        o = o_min;
    else
        o = o_min + (o_max - o_min) * (freqset - freqset(1))./(freqset(end) - freqset(1));
    end
    
    o_int = floor(o);
    al = o - o_int;
    lnC = zeros(length(freqset), length(x));
    
    for i=1:o_max
        idx = i <= o_int;
        freqset_i = freqset(idx);
        
        C_i = cwt_mmorlet(x, fs, w0*i, freqset_i, gpumode);
        
        lnC(idx, :) = lnC(idx, :) + log(C_i);
        
        idx = o_int == i;
        
        if ~all(idx == 0)
            freqset_j = freqset(idx);
            C_i = cwt_mmorlet(x, fs, w0*(i + 1), freqset_j, gpumode);

            al_j = al(idx);
            lnC(idx, :) = lnC(idx, :) + bsxfun(@times, al_j, log(C_i));
        end
    end
    
    C = exp(bsxfun(@times, 1./o, lnC));
    
    P = abs(C).^2;
    
    %%
    figure(2);
    for k=1:K
        subplot(K, 1, k);
        plot(t, 2.*sqrt(P(k, :))); hold on;
        plot([t(1) t(end)], [A(k) A(k)], '-.m'); hold off;
        title(sprintf('f = %3.3f', f(k)));
        ylim([A(k)-0.15 A(k)+0.15]);
    end
end