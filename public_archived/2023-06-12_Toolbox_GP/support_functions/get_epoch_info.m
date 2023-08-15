function epoch = get_epoch_info(ts,i_ep,CE,P,C,for_plot)
%% GET_EPOCH_INFO: returns epoch phase information and filtered data
% Inputs:
%   - ts    - epoch timeseries (e.g. X(channel)(:,i_ep))
%   - i_ep  - index of epoch in in X, for finding epoch events
%   - CE    - channel event struct (SE.(channel).events.Events)
%   - P     - parameter struct from Spectral Events batch
%   - C     - parameter struct from BEAPP data conversion
%   - for_plot - logical, whether to rescale data for visualizing
% Outputs:
%   - epoch - struct with epoch information, including fields:
%       - raw
%       - cwt
%       - phase
%       - filt
%       - event_range
%%

% get optional inputs
if ~exist('for_plot','var'), for_plot = 0; end

% reshape timeseries
epoch.raw = reshape(ts,[1 length(ts)]);
% rescale to range for plotting
if for_plot
    epoch.raw = normalize(epoch.raw,'range',[2 4]);
end
% apply cwt to get phase and filtered series
[epoch.cwt,epoch.cwt_f] = cwt( epoch.raw,'amor',C.Fs, ...
                               FrequencyLimits=P.band.range );
epoch.phase = angle(epoch.cwt);
% get filtered epoch waveform
epoch.filt  = real(epoch.cwt);
% find epoch specific events
epoch.event_range = find([CE.trialind]==i_ep);

end

%% helper functions
function ep = find_epoch_extrema(ep,C)  
% find troughs indeces and their prominence in filtered data
[ep.tr.i,ep.tr.p] = islocalmin(ep.filt,'MinSeparation',C.min_pd);
% find peaks and their prominence in filtered data
[ep.pk.i,ep.pk.p] = islocalmax(ep.filt,'MinSeparation',C.min_pd);
% combine troughs and peaks into extrema
ep.extrema.i = ep.tr.i | ep.pk.i;
ep.extrema.x = C.t_win.samp.i(ep.extrema.i);
ep.extrema.p = ep.tr.p + ep.pk.p;
end