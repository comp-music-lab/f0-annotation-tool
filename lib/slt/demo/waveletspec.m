function maxwid = waveletspec(conf, createfig)
    %% read configuration
    name = conf.name;
    freqset = conf.freqset;
    w0 = conf.w0;
    o_min = conf.o_min;
    o_max = conf.o_max;
    
    %% create impulse input
    fs = 22050;
    T = 6;
    x = zeros(T*fs, 1);
    x(1) = 1;
    
    %%
    ome = (2*pi) .* freqset;
    df = (2*pi) * (fs/length(x));
    F = (0:(length(x) - 1)) .* df;
    
    o = adaptiveorder(freqset, o_min, o_max);
    o = ceil(o);
    
    %%
    dt = 1/fs;
    t = ((0:(length(x) - 1)) - ceil(length(x)/2)) .* dt;
    
    %%
    t_wid = zeros(length(freqset), 1);
    
    for i=1:length(ome)
        H = mmorlet_f(ome(i), w0*o(i), F);
        h_i = ifft(H) ./ dt;
        
        idx = ceil(length(x)/2);
        h_i = [h_i((idx + 1):end) h_i(1:idx)];
        
        h = mmorlet_t(ome(i), w0*o(i), t);
        assert(all(abs(real(h) - real(h_i)) < 1e-6) && all(abs(imag(h) - imag(h_i)) < 1e-6), 'Check wavelet construction');
        
        E = sum(abs(h_i).^2);
        cdf = cumsum(abs(h_i).^2)./E;
        t_left = find(cdf <= 0.005, 1, 'last');
        t_right = find(cdf >= 0.995, 1, 'first');
        
        t_wid(i) = (t(t_right) - t(t_left))/2;
    end
    
    maxwid = max(t_wid);
    
    %%
    if createfig
        figure(1);

        subplot(2, 1, 1);
        plot(freqset, t_wid);
        title(['\omega_{0} = ' num2str(w0), ' (type = ', name, ')']);

        subplot(2, 1, 2);
        plot(freqset, o);
    end
end