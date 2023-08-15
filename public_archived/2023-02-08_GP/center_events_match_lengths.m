function evt = center_events_match_lengths(evt)
    drives = {'event','distal','proximal'};
    n_evt  = length(evt);
    %% center events via zero padding at beginning
    max_pk_i = max([evt(:).dist_pk_i]);
    for i_e = 1:n_evt
        ph_lag = max_pk_i - evt(i_e).dist_pk_i;
        for i_d = 1:length(drives)
            dr = drives{i_d};
            if ~isfield(evt,dr), continue; end
            % add zero padding st current distal peak matches distal peak 
            % with highest latency
            evt(i_e).(dr) = zero_pad(evt(i_e).(dr),ph_lag,0);
        end
        evt(i_e).dist_pk_i = max_pk_i; 
    end

    %% match lengths via zero padding at end
    max_len = max([evt(:).length]);
    for i_e = 1:n_evt
        len_diff = max_len - length(evt(i_e).event);
        for i_d = 1:length(drives)
            dr = drives{i_d};
            if ~isfield(evt,dr), continue; end
            % add zero padding st current length matches longest event
            evt(i_e).(dr) = zero_pad(evt(i_e).(dr),0,len_diff);
        end
        evt(i_e).length = length(evt(i_e).event); 
    end
end