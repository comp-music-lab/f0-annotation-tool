function [o, o_int, al] = adaptiveorder(freqset, o_min, o_max)
    o = o_min + (o_max - o_min) * (freqset - freqset(1))./(freqset(end) - freqset(1));
    
    o_int = floor(o);
    
    al = o - o_int;
end