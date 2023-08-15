function idxs = get_segment_times(seg_win, srate)
    % get times of samples for a given segment window and samp rate
    n_samps = ceil(srate*sumabs(seg_win)/1000);
    idxs = linspace(seg_win(1),seg_win(2),n_samps);
end