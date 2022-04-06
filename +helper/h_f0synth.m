function f0synth = h_f0synth(f0vec, fs, dec_n)
    %%    
    N = round(length(f0vec) + (dec_n - 1)*(length(f0vec) - 1));
    f0vec = interpft(f0vec, N);
    
    %%
    flowTime = cumsum(f0vec./fs);
    f0synth = 0.5 .* sin(2*pi.*flowTime);
end