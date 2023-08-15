function evts = segment_beta_events(data,Fs,plots)
%% SEGMENT_BETA_EVENTS: segments series of beta events into individuals
%   Inputs:
%       - data  - double array | continuous timeseries of beta events
%       - Fs    - (optional) int | sampling rate in Hz, default: 1000
%       - plots - (optional) logical | 1 = plot outputs, 0 = no plots
%   Outputs:
%       - evts  - struct | segmented beta events, and separated laminar
%                 drives
%%
% verify inputs
if ~exist('Fs','var'), Fs = 1000; end
if ~exist('plots','var'), plots = 0; end

%% find extrema
min_pk_dist = 33 * Fs/1000;
trs  = find(islocalmin(data,'MinSeparation',min_pk_dist)); 
pks  = find(islocalmax(data,'MinSeparation',min_pk_dist));
n_tr = length(trs);

%% segment individual events
% create event counter and track last event end index
i_t = 1; evt_i = 1; evt_strt = 1; % track last_end instead if no overlap
while i_t <= n_tr
    %% find distal drive boundaries: surrounding peaks
    tr = trs(i_t);
    last_pk = max(pks(pks < tr));
    next_pk = min(pks(pks > tr));
    % if no next peak, reached end of events
    if isempty(next_pk), break; end
    % if trough not preceded by peak, skip trough
    if isempty(last_pk)
        i_t = i_t + 1; % last_end = floor( (tr+next_pk)/2 ); 
        continue
    end

    %% find event end: after next peak
    next_tr = min(trs(trs > next_pk));
    if ~isempty(next_tr)
        % if another trough, event ends between next peak and trough
        evt_end = floor( (next_pk + next_tr)/2 );
    elseif data(end) > 0
        % if end of data is positive, end event there
        evt_end = length(data);
    else
        % otherwise, reached end of events
        break
    end
    if i_t > 1
        % if not first trough, event starts bw end last trough and peak
        last_tr = max(trs(trs < last_pk));
        if isempty(last_tr), i_t = i_t + 1; continue; end
        evt_strt = floor( (last_tr+last_pk)/2 ); 
    else
        % else, event starts at start of data
        evt_strt = 1;
    end
    evt = normalize(data(evt_strt:evt_end),'center'); % 0-center data

    %% separate laminar drives
    % distal: between last and next peaks
    d_strt = last_pk-evt_strt+1;
    dist = evt(last_pk-evt_strt+1:d_strt+next_pk-last_pk); %data(last_pk:next_pk);
    % zero-pad to have equal array lengths
    dist = zero_pad(dist,last_pk-evt_strt,evt_end-next_pk);
    dist(dist~=0) = dist(dist~=0) - max(dist); % bl normalize
    % proximal: event minus distal
    prox = evt - dist; 
    
    %% save data to structs and update counters
    evts(evt_i).event     = evt;
    evts(evt_i).distal    = dist;
    evts(evt_i).proximal  = prox;
    % track distance to distal trough & evt length for centering later
    evts(evt_i).dist_pk_i = tr-evt_strt+1;
    evts(evt_i).length    = length(evt);
    % save event parameters
    evts(evt_i).distal_amp    = min(dist);
    evts(evt_i).distal_w_ms   = (next_pk - last_pk) * 1000/Fs;
    evts(evt_i).proximal_amp  = max(prox);
    evts(evt_i).proximal_w_ms = length(evt) * 1000/Fs;
    % update counters (if do not want overlap, i_t = i_t + 2)
    i_t = i_t + 1; evt_i = evt_i + 1; last_end = evt_end;
end
if ~exist('evts','var'), evts = []; return; end

%% center all data around distal peak & equate lengths
evts = center_events_match_lengths(evts);

%% plot outputs
if plots
    figure
    subplot(1,n_dr+1,1), plot(data); title('data')
    drives = {'event','distal','proximal'};
    for i_d = 1:length(drives)
        drive = drives{i_d};
        subplot(1,n_dr+1,i_d+1); hold on
        title(drive)
        for i_e = 1:n_evt
            plot(evts(i_e).(drive),'Color','#ccc')
        end
        mean_evt = mean(padcat(evts(:).(drive)),1,'omitnan');
        plot(mean_evt,'LineWidth',2,'Color','#000')
    end
end

end