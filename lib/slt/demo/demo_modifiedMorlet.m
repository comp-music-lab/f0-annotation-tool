function demo_modifiedMorlet
    %%
    dt = 0.001;
    fs = 1/dt;
    T = 8 + round(rand, 3);
    t = -T:dt:T;
    
     %%
    A = rand;
    f = rand * 80;
    x = A.*sin(2*pi*f.*t);
    E = sum(x.^2)/length(x);
    
    %%
    ome = 2 * pi * f;
    c = max(4, rand * 16);
    
    hfun = @(ome, c, t) 5*ome/(c * (2*pi)^1.5) .* exp(-0.5.*(5.*ome.*t/(2*pi*c)).^2) .* exp(1i.*ome.*t);
    Hfun = @(ome, c, f) exp(-(sqrt(2)*pi*c*(f - ome)./(5*ome)).^2);
    C = c*2*sqrt(pi)/(5*ome);
    
    %%
    h = hfun(ome, c, t);
    
    F = (0:(length(t) - 1)) .* 1/(dt*length(t));
    H = Hfun(ome, c, 2*pi.*F);
    
    %%
    idx = find(t == 0);
    h_i = [h(idx:end) h(1:(idx - 1))];
    h_j = ifft(H) ./ dt;
    
    H_i = fft(h_i) .* dt;
    
   %%
    L1norm = sum(abs(h)) * dt;
    
    df = 2*pi * (fs/length(H)); % angular ferquency version
    AC = sum(abs(H).^2) * C * df;
    
    %%
    C = ifft(fft(x) .* H);
    P = abs(C).^2;
    
    %%
    figure(1);
    
    subplot(4, 1, 1);
    plot(t, real(h)); hold on;
    plot(t, imag(h), '-.m'); hold off;
    title(sprintf('c = %3.3f (L1 norm = %e)', c, L1norm));
    xlim([t(find(real(h) > 1e-3, 1, 'first')) t(find(real(h) > 1e-3, 1, 'last'))]);
    
    subplot(4, 1, 2);
    plot(t, real(h_i) - real(h_j)); hold on;
    plot(t, imag(h_i) - imag(h_j), '-.m'); hold off;
    
    subplot(4, 1, 3);
    [~, idx] = max(abs(H).^2);
    plot(F, abs(H).^2);
    title(sprintf('f = %3.3f (peak = %3.3f, AC = %3.3f)', f, F(idx), AC));
    
    subplot(4, 1, 4);
    plot(F, real(H) - real(H_i)); hold on;
    plot(F, imag(H) - imag(H_i), '-.m'); hold off;
    
    figure(2);
    
    subplot(3, 1, 1);
    plot(t, P);
    xlim([t(1) t(end)])
    
    subplot(3, 1, 2);
    plot(t, 2.*P); hold on;
    plot([t(1) t(end)], [E E], '-.m'); hold off;
    title(sprintf('E = %3.3f', E));
    xlim([t(1) t(end)])
    
    subplot(3, 1, 3);
    plot(t, 2.*sqrt(P)); hold on;
    plot([t(1) t(end)], [A A], '-.m'); hold off;
    title(sprintf('A = %3.3f', A));
    xlim([t(1) t(end)])
end

%{
f_i = rand;
x = exp(-(sqrt(2)*pi*c/(5*ome))^2 * (f_i - ome)^2);

sgm = 5*ome/(2*pi*c);
y = exp(-(f_i - ome)^2/(2*sgm^2));

z = normpdf(f_i, ome, sgm) * sqrt(2*pi*sgm^2);

q = Hfun(ome, c, f_i);

disp([x y z q]);
%}