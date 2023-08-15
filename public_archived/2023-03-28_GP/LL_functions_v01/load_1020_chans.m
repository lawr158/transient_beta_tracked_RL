function [chans_1020,chan_layout] = load_1020_chans()
    chans_1020 = [22 9 33 24 124 122 36 104 58 62 96 70 83 45 108 52 92 11];
    net = load_hcgsn129(); chan_layout = net(chans_1020);
end