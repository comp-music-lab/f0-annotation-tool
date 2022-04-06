function demo_cepsmap
    fs = 6000;
    dt = 1/fs;
    N = 8192;
    
    x_F = (0:(N - 1)) .* fs/N;
    x_T = [Inf N*dt./(1:(N - 1))];
    
    figure(1);
    subplot(2, 1, 1);
    plot(x_T);
    title('Time domain (period)');
    set(gca, 'YScale', 'log');
    subplot(2, 1, 2);
    plot(x_F);
    title('Frequency domain');
end