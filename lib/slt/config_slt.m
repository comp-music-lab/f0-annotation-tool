function sltconf = config_slt
    freqset = 440 .* 2.^((-530:500)./120);
    
    name = 'superlet';
    %{
    w0 = 4;
    o_min = 8;
    o_max = 32;
    order = 'f';
    %}
    %{
    w0 = 32;
    o_min = 1;
    o_max = 4;
    order = 'f';
    %}
    %{
    w0 = 16;
    o_min = 12;
    o_max = 24;
    order = 'i';
    %}
    %%{
    w0 = 8;
    o_min = 1;
    o_max = 8;
    order = 'f';
    %}

    N_max = 131072 * 1024;
    
    sltconf = struct('name', name, 'freqset', freqset, 'w0', w0, 'o_min', o_min, 'o_max', o_max, 'N_max', N_max, 'order', order);
end