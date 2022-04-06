function [o, o_int, al] = iadaptiveorder(freqset, o_min, o_max)
    if numel(freqset) == 1
        o = o_min;
    else
        o = o_min + (o_max - o_min) * (log(freqset) - log(freqset(1)))./(log(freqset(end)) - log(freqset(1)));
    end

    o = flipud(o(:));

    o_int = floor(o);
    
    al = o - o_int;
end