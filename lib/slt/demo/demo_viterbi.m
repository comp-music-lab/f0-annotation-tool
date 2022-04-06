function demo_viterbi
    J = 9;
    L = 8;
    d = 2;
    A = rand(J, L);
    
    f_start = randi(J);
    f_end = randi(J);
    
    y = algorithm1(A, f_start, f_end, d);
    cost_y = pathcost(A, y, d, f_start);
    
    z = algorithm2(A, f_start, f_end, d);
    assert(all(all(z == viterbidec(A, f_start, f_end, d))), 'check viterbidec')
    cost_z = pathcost(A, z, d, f_start);
    
    %%{
    tic;
    q = fullsearch(A, f_start, f_end, d);
    cost_q = pathcost(A, q, d, f_start);
    toc;
    %}
    
    figure(1);
    subplot(4, 1, 1);
    imagesc(1:L, 1:J, A);
    view(0, 90);
    axis tight;
    
    subplot(4, 1, 2);
    imagesc(1:L, 1:J, y);
    view(0, 90);
    axis tight;
    title(sprintf('cost = %3.3f', cost_y));
    
    subplot(4, 1, 3);
    imagesc(1:L, 1:J, z);
    view(0, 90);
    axis tight;
    title(sprintf('cost = %3.3f', cost_z));
    
    %%{
    subplot(4, 1, 4);
    imagesc(1:L, 1:J, q);
    view(0, 90);
    axis tight;
    title(sprintf('cost = %3.3f', cost_q));
    %}
end

function y = fullsearch(x, f_start, f_end, d)
    z = x .* 0;
    
    cost = Inf;

    J = size(x, 1);
    L = size(x, 2);
    
    z(f_start, 1) = x(f_start, 1);
    z(f_end, L) = x(f_end, L);

    counter = ones(1, L - 2);
    
    while ~all(counter == J)
        z(:, 2:(L - 1)) = 0;
        
        for i=2:(L - 1)
            z(counter(i - 1), i) = x(counter(i - 1), i);
        end

        cost_t = pathcost(x, z, d, f_start);
        
        if cost_t < cost
            cost = cost_t;
            y = z;
        end
        
        for i = (L - 2):-1:1
            counter(i) = counter(i) + 1;
            
            if counter(i) > J
                counter(i) = 1;
            else
                break;
            end
        end
    end
end

function cost = pathcost(x, y, d, f_start)
    cost = 0;
    J = size(x, 1);
    L = size(x, 2);
    c = 2.5;
    
    m = f_start;
    
    for i=2:L
        l = find(y(:, i) ~= 0);
        
        g = abs(l - m);
        if g <= d
            g = 0;
        else
            g = c*(g - d);
        end
        
        f = J - tiedrank(x(:, i));
        
        cost = cost + f(l) + g;
        
        m = l;
    end
end

function y = algorithm1(x, f_start, f_end, d)
    J = size(x, 1);
    L = size(x, 2);
    y = x .* 0;
    
    c = 2.5;
    I = (1:J)';
    
    y(f_end, L) = x(f_end, L);
    y(f_start, 1) = x(f_start, 1);
    k = f_end;
    
    for i=(L-1):-1:2
        f = J - tiedrank(x(:, i));
        
        g = abs(I - k);
        g(g <= d) = 0;
        g(g > d) = c.*(g(g > d) - d);
        
        [~, k_i] = min(f + g);
        
        y(k_i, i) = x(k_i, i);
        k = k_i;
    end
end

function y = algorithm2(x, f_start, f_end, d)
    J = size(x, 1);
    L = size(x, 2);
    y = x .* 0;
    
    c = 2.5;
    I = (1:J)';
    
    %%
    M = zeros((L - 1)*J, 5);
    
    %%
    M(1:J, 1:3) = repmat([1 f_start 2], J, 1);
    M(1:J, 4) = 1:J;
    
    g = abs(I - f_start);
    g(g <= d) = 0;
    g(g > d) = c.*(g(g > d) - d);
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

                g = abs(I - j);
                g(g <= d) = 0;
                g(g > d) = c.*(g(g > d) - d);

                cost = g + p_tmp;
                cost(1:(idx_min - 1)) = Inf;
                cost((idx_max + 1):end) = Inf;

                [~, l_tmp] = min(cost);

                pi_l = M(M(:, 1) == (i - 1) & M(:, 4) == l_tmp, 5);
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