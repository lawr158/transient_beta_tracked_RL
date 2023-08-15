function [rel_pk,w,wconv] = anaylze_beta_event(event,Fs,plots)
%% ANALYZE_BETA_EVENT: analyzes waveform input beta event
%   Inputs:
%       - event - double array | time series with beta event data
%       - Fs    - int | sampling rate of data in Hz, default: 1000
%       - pltos - logical | 1 = generate plots, 0 = no plots, default: 0
%   Outputs:
%       - rel_pk - struct | peaks detected by each wavelet, in relative amp
%       - w      - struct | wavelets used for each laminar drive
%       - wconv  - struct | results of convolution with each wavelet 
%% set up for analysis
% verify inputs
if ~exist('Fs','var'), Fs = 1000; end
if ~exist('plots','var'), plots = 0; end

% get event duration and generate beta event wavelets
evt_dur = length(event) * 1000/Fs;
d.dur   = evt_dur/2; p.dur = evt_dur;
[~,x,w.d,w.p,w.b] = sim_beta_event(d,p,Fs,[]);

%% analyze event
% convolve data with separated laminar drives
evt_half = round(length(event)/2);
wconv.d = conv(event,w.d); wconv.d = wconv.d(evt_half:end-evt_half);
wconv.p = conv(event,w.p); wconv.p = wconv.p(evt_half:end-evt_half);
wconv.b = conv(event,w.b); wconv.b = wconv.b(evt_half:end-evt_half);

% find loxal maxima of distal/proximal convolutions
peaks_d = islocalmax(wconv.d); 
peaks_p = islocalmax(wconv.p);
peaks_b = islocalmax(wconv.b); peaks_b = peaks_b - islocalmin(wconv.b);

% find relative amplitudes of maxima
rel_amp  = rescale([wconv.d;wconv.p;wconv.b],-1,1);
rel_pk.d = peaks_d .* rel_amp(1,:);
rel_pk.p = peaks_p .* rel_amp(2,:);
rel_pk.b = peaks_b .* rel_amp(3,:);

%% plot convolution results & maxima
if plots
    % plot event(s) % laminar drive wavelets
    figure; hold on
    plot(x,event,'LineWidth',2,'Color','#000'); 
    plot(x,w.d); plot(x,w.p); plot(x,w.b); 
    legend({'event','distal drive','proximal drive','both drive'})
    title('Beta event & wavelets')
    xlabel('Time (ms)'); ylabel('Amplitude') 
    
    % get x axis values for convolution
    wc_len = length(wconv.d);
    x_wc   = linspace(0,length(wconv.d)*1000/Fs,wc_len);
    t_wc   = 1:wc_len; t_ss = 1:ceil(wc_len/10):wc_len;

    % plot convolution timeseries
    figure; hold on
    plot(x_wc,wconv.d); plot(x_wc,wconv.p); plot(x_wc,wconv.b);
    legend({'distal','proximal','both'})
    title('Convolution: event x laminar drives'); 
    xlabel('Time (ms)'); ylabel('Amplitude')

    % plot convolution heatmap
    figure; hold on
    pcolor([wconv.d;wconv.p;wconv.b]); colorbar; shading interp
    title('Convolution: event x laminar drives'); 
    xlabel('Time (ms)'); ylabel('Wavelet');
    xticks(t_wc(t_ss)); xticklabels(round(x_wc(t_ss)))
    yticks([1:3]); yticklabels({'Distal','Proximal','Both'})

    % plot convolution maxima
    figure; hold on
    plot(x_wc,rel_pk.d); plot(x_wc,rel_pk.p); plot(x_wc,rel_pk.b)
    title('Detected Peaks'); legend({'distal','proximal','both'})
    xlabel('Time (ms)'); ylabel('Relative amplitude')
end

end