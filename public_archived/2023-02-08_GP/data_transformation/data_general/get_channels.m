function channels = get_channels(brain_regs)
    if isempty(brain_regs)
        channels = [];
        return
    end
    
    channmap = containers.Map({'C3','C4','F3','F4','F7','F8','Fp1',...
        'Fp2','Fz','O1','O2','Pz','P3','P4','P7','P8','T5','T6',...
        'T3','T4','T7','T8'},...
        {36,104,24,124,33,122,22,9,11,70,83,62,52,92,58,96,58,96,...
        45,108,45,108});
    
    n_chans  = size(brain_regs,2);
    channels = zeros(1,n_chans);
    for i=1:n_chans
        channels(i) = channmap(brain_regs{i});
    end
end