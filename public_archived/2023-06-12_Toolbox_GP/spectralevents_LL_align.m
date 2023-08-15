function SE_ts = spectralevents_LL_align(SE_tfr,X,P,C)
%% SPECTRALEVENTS_LL_ALIGN: align TFR events to TS in channel data
% Inputs:
%   - SE_tfr    - TFR spectralevents output struct for a specific channel
%   - X         - timeseries data (in spectralevents format) for the same
%                 channel as SE_tfr
%   - P         - spectralevents batch parameters struct
%   - C         - BEAPP to spectralevents conversion parameters struct
% Outputs:
%   - SE_ts - struct with aligned event information. aligned data are added
%             on top of SE_tfr.events.Events
%
% find range of cycle periods for extrema finding/event centering
pd_rng = 1000./P.band.range * C.Fs/1000;
% get channel events struct
SE_ts = SE_tfr.events.Events;
% iterate through epochs
n_epoch = size(X.data,2);
for i_ep = 1:n_epoch
    % get epoch filtered ts, extrema, event range
    epoch = get_epoch_info(X.data(:,i_ep),i_ep,SE_ts,P,C); 
    % iterate through events
    for i_ev = epoch.event_range
        % get ts event information
        evt = align_tfr_event_ts(SE_ts(i_ev),epoch,C.Fs,C.t_win,pd_rng);
        % add ts fields to channel struct
        ts_fields = fieldnames(evt);
        for i_f = 1:length(ts_fields)
            tsf = ts_fields{i_f};
            SE_ts(i_ev).(tsf) = evt.(tsf);
        end
    end
end

end