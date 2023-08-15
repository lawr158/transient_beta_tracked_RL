function chan_str = get_chan_str(net_1020_chans, other_chans)
    net1020_str = strjoin(net_1020_chans, '');
    
    other_str = '';
    for i = 1:size(other_chans,2)
        other_str = [other_str 'E' int2str(other_chans(i))];
    end

    chan_str = [net1020_str other_str];
end