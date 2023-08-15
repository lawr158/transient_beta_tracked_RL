function [chan_nos,chan_layout] = load_1020_chans(net)

% check input, set default if none
if ~exist('net','var'), net = 'hcgsn129'; end

switch net
    case 'hcgsn129'
    chan_nos = [22 9 33 24 124 122 36 104 58 62 96 70 83 45 108 52 92 11];
    [~,chan_layout] = load_hcgsn129;
end

end