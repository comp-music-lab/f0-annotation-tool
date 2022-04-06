function demo_fslt3
    %%
    %{
    audiofile = '../audio/AmaSua.m4a';
    t_start = 4;
    t_end = 20;
    %}
    
    %{
    audiofile = '../audio/CantoMetConsensus_02.mp3';
    t_start = 4;
    t_end = 19.8;
    %}
    
    %{
    audiofile = '../audio/CantoMetConsensus_05.mp3';
    t_start = 4;
    t_end = 18;
    %}
    
    %{
    audiofile = '../audio/Hirado-bushi.wav';
    t_start = 15.2;
    t_end = 28;
    %}
    
    %{
    audiofile = '../audio/PianoSonataNo27_Beethoven.m4a';
    t_start = 0;
    t_end = 9;
    %}
    
    %{
    audiofile = '../audio/Fandango_Rodrigo.m4a';
    t_start = 0.5;
    t_end = 8;
    %}
    
    %{
    audiofile = '../audio/CantoMetConsensus_17.mp3';
    t_start = 8.7;
    t_end = 19.9;
    %}
    
    %%{
    audiofile = '../audio/NAIV-075.wav';
    t_start = 0;
    t_end = 14;
    %}
    
    %{
    audiofile = '../audio/NHSDiscography-050.wav';
    t_start = 6.5;
    t_end = 14.1;
    %}
    
    %{
    audiofile = '../audio/NHSDiscography-092.wav';
    t_start = 2.5;
    t_end = 14.5;
    %}
    
    %{
    audiofile = '../audio/OsJusti_Bruckner.m4a';
    t_start = 0;
    t_end = 20;
    %}
    
    %%
    figsize = round([1 (1+sqrt(5))/2] .* 500);
    
    fs_thresh = 16000;
    dec_t = 0.005;
    
    gpumode = true;
    freqset = 440 .* 2.^((-530:500)./120);
    
    gam_X = 0.01;
    gam_Y = 1.0;
    shrapness = 1;

    %{
    w0 = 8;
    o_min = 1;
    o_max = 32;
    prefix = '_superlet';
    %}
    
    %%{
    w0 = 24;
    o_min = 1;
    o_max = 1;
    prefix = '_wavelet';
    %}
    
    %% Read audio
    [x, fs] = audioread(audiofile);
    
    if size(x, 2)
        x = mean(x, 2);
    end
    
    x = x(:)';
    
    %% Segmentation
    n_start = round(fs*t_start + 1);
    n_end = round(fs*t_end);
    x = x(n_start:n_end);
    
    t = t_start + (0:(length(x) - 1)) ./ fs;
    
    %% Downsampling
    while (fs/2) > fs_thresh
        fs = fs/2;
        x = x(1:2:end);
        t = t(1:2:end);
    end
    
    freqset(freqset > fs/2) = [];
    
    %% Superlet
    o = o_min + (o_max - o_min) * (freqset - freqset(1))./(freqset(end) - freqset(1));
    o_int = floor(o);
    al = o - o_int;
    
    dec_n = round(dec_t * fs);
    lnC = zeros(length(freqset), length(x(1:dec_n:end)));
    
    if gpumode
        lnC = gpuArray(cast(lnC, 'single'));
        al = gpuArray(cast(al, 'single'));
        o = gpuArray(cast(o, 'single'));
    end
    
    for i=1:o_max
        idx = i <= o_int;
        freqset_i = freqset(idx);
        
        C_i = cwt_mmorlet(x, fs, w0*i, freqset_i, gpumode);
        
        lnC(idx, :) = lnC(idx, :) + log(C_i(:, 1:dec_n:end));
        
        idx = o_int == i;
        
        if ~all(idx == 0)
            freqset_j = freqset(idx);
            C_j = cwt_mmorlet(x, fs, w0*(i + 1), freqset_j, gpumode);
            C_j = C_j(:, 1:dec_n:end);

            al_j = al(idx)';
            lnC(idx, :) = lnC(idx, :) + bsxfun(@times, al_j, log(C_j));
        end
    end
    
    C = exp(bsxfun(@times, 1./o', lnC));

    P = abs(C).^2;
    W = 2.*sqrt(P);
    
    t_dec = t(1:dec_n:end);
    
    %% Gaussian decompositioin
    %{
    thresh = 0.01;
    d = 6;
    
    buf = gather(2.*P);
    R = buf .* 0;
    
    for i=1:size(R, 2)
        X = buf(:, i);
        
        Y = gaussclean(X, thresh, d);
        R(:, i) = Y;
    end
    %}
    
    %% CFP representation
    CFP = cfprep_log(2.*P, freqset, gam_X, gam_Y, shrapness);
    
    %%
    f = figure(1);
    
    surf(t_dec, freqset, W, 'edgecolor', 'none');
    view(0, 90);
    axis tight;
    ylim([-Inf freqset(end)]);
    set(gca, 'YScale', 'log');
    
    f.Position = [1 1 figsize(2) figsize(1)];
    
    f = figure(2);
    
    surf(t_dec, freqset, log(2.*P), 'edgecolor', 'none');
    view(0, 90);
    axis tight;
    ylim([-Inf freqset(end)]);
    set(gca, 'YScale', 'log');
    
    f.Position = [1 1 figsize(2) figsize(1)];
    
    %{
    f = figure(3);
    
    surf(t, freqset, log(R + min(R(R ~= 0))), 'edgecolor', 'none');
    view(0, 90);
    axis tight;
    ylim([-Inf freqset(end)]);
    set(gca, 'YScale', 'log');
    
    f.Position = [1 1 figsize(2) figsize(1)];
    %}
    
    f = figure(4);
    
    surf(t_dec, freqset, log(CFP + min(CFP(CFP ~= 0))), 'edgecolor', 'none');
    view(0, 90);
    axis tight;
    ylim([-Inf freqset(end)]);
    set(gca, 'YScale', 'log');
    
    f.Position = [1 1 figsize(2) figsize(1)];
    
    %%
    V = W .* 0;
    
    fprintf('Click! (start -> end -> top -> bottom)\n');
    figure(1);
    [xpos, ypos] = ginput(4);
    
    [~, t_lb] = min(abs(t_dec - xpos(1)));
    [~, t_ub] = min(abs(t_dec - xpos(2)));
    [~, f_start] = min(abs(freqset - ypos(1)));
    [~, f_end] = min(abs(freqset - ypos(2)));
    [~, f_ub] = min(abs(freqset - ypos(3)));
    [~, f_lb] = min(abs(freqset - ypos(4)));
    
    V_i = viterbidec(gather(W(f_lb:f_ub, t_lb:t_ub)), f_start - f_lb + 1, f_end - f_lb + 1, 6);
    V(f_lb:f_ub, t_lb:t_ub) = V_i;
   
    %%
    f0vec = zeros(1, size(V, 2));
    for i=1:size(V, 2)
        idx = find(V(:, i) ~= 0);
        
        if ~isempty(idx)
            f0vec(i) = freqset(idx);
        end
    end
    
    f0vec_tmp = f0vec(t_lb:t_ub);
    N_tmp = length(f0vec_tmp) + (dec_n - 1)*(length(f0vec_tmp) - 1);
    f0vec_tmp = interpft(f0vec_tmp, N_tmp);
    
    flowTime = cumsum(f0vec_tmp./fs);
    f0synth = 0.5 .* sin(2*pi.*flowTime);
    
    [~, tt_start] = min(abs(t - t_dec(t_lb)));
    [~, tt_end] = min(abs(t - t_dec(t_ub)));
    x_tmp = x(tt_start:tt_end);
    
    sound(f0synth, fs);
    sound(x_tmp, fs);
    sound(f0synth + x_tmp, fs);
    
    %%
    s = strsplit(audiofile, '/');
    s = strsplit(s{end}, '.');
    
    figfilename = strcat('./figure/', s{1}, prefix, '_amp.png');
    saveas(figure(1), figfilename);
    
    figfilename = strcat('./figure/', s{1}, prefix, '_pow.png');
    saveas(figure(2), figfilename);
    
    %{
    figfilename = strcat('./figure/', s{1}, prefix, '_peak.png');
    saveas(figure(3), figfilename);
    %}
    
    figfilename = strcat('./figure/', s{1}, prefix, '_cfp.png');
    saveas(figure(4), figfilename);
end