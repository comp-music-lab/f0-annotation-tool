function H = mmorlet_f(ome, c, f)
    H = exp(-(sqrt(2)*pi*c*(f - ome)./(5*ome)).^2);
end