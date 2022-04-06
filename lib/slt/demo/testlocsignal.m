function [x, t, f0, a, A] = testlocsignal(fs, T, K)
    t = (0:(fs*T - 1))./fs;
    x = zeros(1, length(t));
    
    f0 = max(30, rand(K, 1) .* (0.8 * fs/2));
    f0 = sort(f0);
    a = max(rand(K, 1), 0.1);
    n_start = randi(length(t), [K 1]);
    len = arrayfun(@(x) randi(length(t) - x), n_start);
    n_end = n_start + len;
    
    A = zeros(K, length(t));
    
    for k=1:K
        x(n_start(k):n_end(k)) = x(n_start(k):n_end(k)) + ...
            a(k).*sin(2*pi*f0(k).*(t(n_start(k):n_end(k)) - rand)) .* hann(len(k) + 1)';
        
        A(k, n_start(k):n_end(k)) = a(k) .* hann(len(k) + 1)';
    end
end