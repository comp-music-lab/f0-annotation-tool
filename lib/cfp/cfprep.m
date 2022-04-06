function Z = cfprep(X, freqset, gam_X, gam_Y, sharpness)
    %%
    if isa(X, 'gpuArray')
        freqset = gpuArray(cast(freqset, 'single'));
    end
    
    %%
    J = length(freqset);
    
    maxpartial = 8;
    
    Y = X .* 0;
    L = abs(X).^gam_X;
    
    %%
    L(isnan(L)) = 0;
    
    %%
    for j=1:J
        f0 = freqset(j);
        phi = 2.*pi.*(freqset - f0)./f0;
        filter = cos(phi);
        
        filter(phi < -pi/2) = 0;
        
        filter(phi > ((maxpartial + 1)*2*pi)) = 0;
        
        filter = sign(filter) .* abs(filter).^sharpness;

        Y(j, :) = filter * L;
    end
    
    %%
    Y = gam_Y.^Y;

    %%
    %Z = Y.^L;
    Z = Y;
end

%{
figure;
surf(1:size(Y, 2), freqset, Y, 'edgecolor', 'none');
view(0, 90);
axis tight;
set(gca, 'YScale', 'log');
%}

%{
figure;
subplot(2, 1, 1);
plot(freqset, filter); hold on;
stem(f0.*(1:(4 + 1)), ones(1, 4 + 1)); hold off;
title(sprintf('f_0 = %3.3f', f0));
subplot(2, 1, 2);
plot(filter);
%}