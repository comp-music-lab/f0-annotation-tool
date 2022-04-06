function y = viterbidec(x, f_start, f_end, d)
    %%
    J = size(x, 1);
    L = size(x, 2);
    y = x .* 0;
    
    c = 100;
    I = (1:J)';
    
    %%
    M = zeros((L - 1)*J, 5);
    
    %%
    M(1:J, 1:3) = repmat([1 f_start 2], J, 1);
    M(1:J, 4) = 1:J;
    
    g = distfun(I, f_start, d, c);
    f = J - tiedrank(x(:, 2));
    
    M(1:J, 5) = (f + g)';
    
    for i=2:(L - 1)
        offset = J*(i - 1);
        M((offset + 1):(offset + J), 1) = i;
        M((offset + 1):(offset + J), 3) = i + 1;
        M((offset + 1):(offset + J), 4) = 1:J;
        
        p_tmp = M(M(:, 1) == (i - 1), 5);
        pi_tmp = min(p_tmp);
        f = J - tiedrank(x(:, i + 1));
        
        for j=1:J
            rho = d;
            
            while true
                idx_min = max(1, j - rho);
                idx_max = min(J, j + rho);

                g = distfun(I, j, d, c);

                cost = g + p_tmp;
                cost(1:(idx_min - 1)) = Inf;
                cost((idx_max + 1):end) = Inf;

                [~, l_tmp] = min(cost);

                pi_l = p_tmp(l_tmp);
                A = pi_l + g(l_tmp);
                B = c*(rho + 1 - d) + pi_tmp;

                if A < B
                    M((offset + j), 2) = l_tmp;
                    M((offset + j), 5) = pi_l + g(l_tmp) + f(j);

                    break;
                else
                    rho = rho + 1;
                end
            end
        end
    end
    
    %%
    k = f_end;
    
    for i=L:-1:2
        y(k, i) = x(k, i);
        
        backward_state = M(M(:, 3) == i, 2);
        
        k = backward_state(k);
    end
    
    y(f_start, 1) = x(f_start, 1);
end

function g = distfun(I, j, d, c)
    g = abs(I - j);
    
    g(g <= d) = 0;
    
    g(g > d) = c.*(g(g > d) - d);
end

%{
figure;
surf(t_dec, freqset, V.^0.01, 'edgecolor', 'none');
view(0, 90);
axis tight;
set(gca, 'YScale', 'log');
%}