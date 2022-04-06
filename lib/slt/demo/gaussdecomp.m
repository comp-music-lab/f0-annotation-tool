function [mu_k, sgm_k, const_k, ene_k] = gaussdecomp(X, thresh, d)
    E = sum(X);
    E_i = E;
    mu_k = [];
    sgm_k = [];
    const_k = [];
    ene_k = [];
    
    L = length(X);
    I = (1:L)';
    
    counter = 1;
    thresh_ignore = 30;
    
    while E_i/E > thresh
        [~, locs] = findpeaks(X, 'SortStr', 'descend', 'Npeaks', 1, 'SortStr', 'descend');
        
        idx = locs;
        h = X(idx);
        
        idx_l = max(1, idx - d):max(1, idx - 1);
        r_l = X(idx_l')./h;
        x_l = abs(idx_l - idx)';
        
        idx_r = min(L, idx + 1):min(L, idx + d);
        r_r = X(idx_r')./h;
        x_r = (idx_r - idx)';
        
        sgm_lr = [x_l./sqrt(2.*log(1./r_l)); x_r./sqrt(2.*log(1./r_r))];
        sgm_lr = sgm_lr(sgm_lr ~= 0);
        sgm = median(sgm_lr);
        
        const = sqrt(2*pi*sgm^2) * h;
        gaussfit = const.*normpdf(I, idx, sgm);
        
        Y = X - gaussfit;
        Y(Y < 0) = 0;
        
        E_i = sum(Y);
        mu_k(end + 1) = idx;
        sgm_k(end + 1) = sgm;
        const_k(end + 1) = const;
        ene_k(end + 1) = sum(X) - sum(Y);
        
        %{
        figure(4);
        plot(1:length(X), X); hold on;
        Xx = (idx - 30):(idx + 30);
        plot(Xx, const.*normpdf(Xx, idx, sgm), '-.m'); hold off;
        ylim([0 const_k(1)/sqrt(2*pi*sgm_k(1)^2)*1.05]);
        drawnow;
        pause(0.5);
        %}
        
        X = Y;
        
        counter = counter + 1;
        if counter >= thresh_ignore
            break;
        end
    end
end