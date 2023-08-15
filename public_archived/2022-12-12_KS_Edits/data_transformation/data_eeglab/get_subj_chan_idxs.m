function chan_idxs = get_subj_chan_idxs(subj_info,chans)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    chan_labels = {subj_info.chanlocs.labels};
    chan_idxs   = [find(contains(chan_labels,chans))];
end