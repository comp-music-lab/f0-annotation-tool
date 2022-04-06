function h = mmorlet_t(ome, c, t)
    h = 5*ome/(c * (2*pi)^1.5) .* exp(-0.5.*(5.*ome.*t/(2*pi*c)).^2) .* exp(1i.*ome.*t);
end