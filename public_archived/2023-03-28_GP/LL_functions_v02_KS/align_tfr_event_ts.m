function [evt,tfr] = align_tfr_event_ts(E,epoch,Fs,t_win,pd_rng)
%% ALIGN_TFR_EVENT_TS: finds temporal info from TFR-found spectral event
% Inputs:
%   - E     - spectral events event struct 
%   - epoch - epoch information struct, from get_epoch_info
%   - Fs    - data sampling rate (Hz)
%   - t_win  - struct with segment information from conversion
%   - pd_rng - range of cycle periods, calculated from analysis freq band
% Outputs:
%   - evt - struct with timeseries event information, including fields:
%       - ts_center_i_event - sample index of evt center (from evt onset)
%       - ts_center_i_epoch - sample index of evt center (from epoch onset)
%       - ts_center_ms - time in ms of evt center (from epoch onset)
%       - ts_onset_i_epoch - sample index of evt onset (from epoch onset)
%       - ts_offset_i_epoch - sample index of evt offset
%       - ts_onset_ms - time in ms of evt onset
%       - ts_offset_ms - time in ms of evt offset
%   - tfr - struct with tfr event information
%       - win - timing of tfr event
%       - i   - sample indeces of event in epoch
%
%% get tfr event timing/epoch sample indices
tfr.win = [E.tfr_onset E.tfr_offset]*1000;
tfr.i = get_time_win_idxs(tfr.win,t_win.to_keep,Fs);
% tfr.y = epoch.filt(tfr.i);
%% get minima overlapping with tfr event
% get frequencies of event to get relvant epoch filter data
freq_i = get_indeces([E.upperboundFspan E.lowerboundFspan],epoch.cwt_f,1);
epoch.filt  = sum(epoch.filt(freq_i,:),1);
epoch.phase = unwrap(epoch.phase(round(median(freq_i)),:));
% find troughs in timeseries at event frequencies
trough = find_troughs(epoch.filt(:,tfr.i),pd_rng(2));
%% find event center: max magnitude trough
[~,center_i_tfr_event] = max(trough.p);
% get event center index in epoch
evt.ts_center_i_epoch = tfr.i(center_i_tfr_event);
evt.ts_center_ms = evt.ts_center_i_epoch/Fs*1000; 
%% get timing of event in time series;
% get phases of event center and bounds (+/- 1.5 periods)
ctr_phase = epoch.phase(evt.ts_center_i_epoch);
phase_bounds = [ctr_phase-3*pi ctr_phase ctr_phase+3*pi];
% get bounds of event in epoch samples and ms 
bounds_i_epoch = get_indeces(phase_bounds,epoch.phase);
evt.ts_center_i_event = bounds_i_epoch(2) - bounds_i_epoch(1);
bounds_t_epoch = bounds_i_epoch/Fs*1000;
evt.ts_onset_i_epoch = bounds_i_epoch(1);
evt.ts_offset_i_epoch = bounds_i_epoch(end);
evt.ts_onset_ms = bounds_t_epoch(1); evt.ts_offset_ms = bounds_t_epoch(2);
evt.ts_duration = evt.ts_offset_ms - evt.ts_onset_ms;
%% get trough-centered time series of event
% initialize zero-padded time series, with max possible event duration in
% frequency band (1.5 cycles on each side of trough)
% - in case first half of event defined by very slow event (and thus is
%   longer than expected max), padding by 2.5 cycles
evt_n_samp = 1 + bounds_i_epoch(end) - bounds_i_epoch(1);
max_evt_half_dur = ceil(2.5 * ceil(pd_rng(1))); 
evt_padding = max_evt_half_dur - evt.ts_center_i_event;
evt.y = nan(1,2*max_evt_half_dur+1); evt.y_filt = evt.y; 
evt_i = [1+evt_padding:evt_padding+evt_n_samp];
% skip event if too long (would happen if second half of event is
% disproportionately long)
if evt_padding + evt_n_samp > 2*max_evt_half_dur+1, return; end
try
    % add event data with padding so that trough is centered
    evt.y(evt_i) = epoch.raw(evt.ts_onset_i_epoch:evt.ts_offset_i_epoch);
    evt.y_filt(evt_i) = epoch.filt(evt.ts_onset_i_epoch:evt.ts_offset_i_epoch);
catch
    disp(['First half of event is much longer than expected: '          ...
          int2str(evt.ts_center_i_event) ' samples']);
end
end

%% helper functions
function trough = find_troughs(ts,min_pd)
% find trough indeces and their prominence in filtered data
[trough.i,trough.p] = islocalmin(ts,'MinSeparation',min_pd);
trough.p(trough.p==0) = NaN;
end