function [SE_ts,SE_tfr,TFRs] = spectralevents_LL( P, C, X )
%% SPECTRALEVENTS_LL Find and analyze transient spectral events in a 
%   time-series dataset. Spectral events are defined as local maxima above 
%   a power threshold of a specified band in the non-averaged 
%   time-frequency responses (TFR). Modified for data from the Levin Lab.
%   Last updated by: Gerardo Parra, 2023-03-29
%
%  [SE,TFRs] = SPECTRALEVENTS(P,C,X)
%   
%   Returns a structure array of spectral event features (specEv_struct),
%   cell array containing the time-frequency responses (TFRs), and cell 
%   array of all time-series trials (X) for each subject/session within the
%   dataset comparing various experimental conditions or outcome states 
%   corresponding to each trial. Note: this function sets the factors of
%   median threshold (thrFOM) at 6.
%
% Inputs:
%   P.band.range - range of frequencies ([Fmin_event Fmax_event]; Hz) over 
%       which above-threshold spectral
%       power events are classified.
%   P.fVec - frequency vector (Hz) over which the time-frequency response 
%       (TFR) is calculated. Note that this set must fall within the range 
%       of resolvable/alias-free frequency values (i.e. Fmin>=1/(trial 
%       duration), Fmax<=(Nyquist freq)).
%   C.Fs - sampling frequency (Hz).
%   P.method - integer value specifying which event-finding method to use
%       (1, 2, or 3). Note that the method specifies how much overlap
%       exists between events. Use 1 to replicate the method used in 
%       et al. eLife 2017.
%   X - struct array with the following fields for each channel i_c
%       - X(i_c).channel - channel name
%       - X(i_c).data - m-by-n matrix representing the time-series trials 
%         of the given channel. m is the number of timepoints and n is the
%         number of trials. Note that m timepoints must be uniform across 
%         all trials and channels.
%       - X(i_c).class - numeric or logical trial classification labels; 
%         associates each trial of the given channel to an experimental 
%         condition/outcome/state (e.g., hit
%   P.vis - logical value that determines whether to run basic feature 
%       analysis and output standard figures. or miss, detect or non-detect, 
%         attend-to or attend away).
%   P.ncyc - optional, integer number of cycles to use for wavelet transform 
%   P.thrFOM - optional, integer to use as FOM threshold
%
% Outputs:
%   SET - array of event feature structures, each corresponding with a 
%         channel, respectively.
%   TFRs - cell array with each cell containing the time-frequency response 
%          (freq-by-time-by-trial) for a given channel.
%
% See also SPECTRALEVENTS_FIND, SPECTRALEVENTS_TS2TFR, SPECTRALEVENTS_VIS.

%   -----------------------------------------------------------------------
%   SpectralEvents::spectralevents
%   Copyright (C) 2018  Ryan Thorpe
%
%   This file is part of the SpectralEvents toolbox.
% 
%   SpectralEvents is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   SpectralEvents is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.
%   -----------------------------------------------------------------------
%% validate inputs
if ~isfield(P,'ncyc'), P.ncyc = 5; end
if ~isfield(P,'thrFOM'), P.thrFOM = 6; end
if ~isfield(P,'feats'), P.feats = { 'ts_duration', 'y', 'y_filt' }; end

% Validate format of trial class labels and reformat if necessary
n_chan = numel(X); % Number of subjects/sessions
for i_c=1:n_chan
    if ~(isnumeric(X(i_c).class) || islogical(X(i_c).class))
        error('Trial classification labels must be numeric or logical.')
    end
    cl_s = X(i_c).class;
    if numel(cl_s)==1
        % Reformat if entered as single values instead of vectors
        cl_s = ones(1,size(X(i_c),2))*cl_s; X(i_c).class = cl_s;
    elseif isequal(size(cl_s),[size(X(i_c).data,2),1])
        X(i_c).class = cl_s';
    elseif ~isequal(size(cl_s),[1,size(X(i_c).data,2)])
        error(['Trial classification labels must be formatted as '      ...
               'either a 1-row array or single value.'])
    end
end
clear cl_s
% Validate P.fVec input
Fn = C.Fs/2; % Nyquist frequency
dt = 1/C.Fs; % Sampling time interval
nsamp = max(cellfun(@(x) size(x,1),{X.data})); % number of samples in data
Fmin = 1/(nsamp*dt); % Minimum resolvable frequency
if P.fVec(1)<Fmin || P.fVec(end)>Fn || abs(P.fVec(2)-P.fVec(1))<Fmin
    error(['Frequency vector includes values outside the resolvable/'   ...
           'alias-free range.'])
end

%% Solve for the time-frequency response (TFR) and run spectral events
%   analysis on each trial within each subject/session
TFRs = struct('channel',{X.channel},'tfr',NaN); 
SE_tfr = struct('channel',{X.channel},'trials',NaN,'events',NaN,'IEI',NaN);
SE_ts  = struct('channel',{X.channel},'events',NaN);
for i_c=1:n_chan
    % skip if nan channel
    if isnan(X(i_c).data), continue; end
    TFR = []; % Matrix for storing freq-by-time-by-trial
    for i_t=1:size(X(i_c).data,2)
        % Transformation calculated using a Morlet wavelet (width=7), see 
        %   4DToolbox for function details)
        [TFR_t,tVec,~] = spectralevents_ts2tfr( X(i_c).data(:,i_t),     ...
                                                P.fVec, C.Fs, P.ncyc ); 
        % Concatenate each trial along the 3rd dimension
        TFR = cat(3,TFR,TFR_t); 
    end
    TFRs(i_c).tfr = TFR; % Append TFR for the given channel
    if any(TFR == 0), continue; end
    % Find spectral events
    [SE_tfr(i_c).trials,SE_tfr(i_c).events,SE_tfr(i_c).IEI] =           ...
        spectralevents_LL_find( P.method, P.band.range, P.thrFOM, tVec, ...
                                P.fVec, TFR, X(i_c).class );
    % align TFR events to timeseries
    SE_ts(i_c).events = spectralevents_LL_align( SE_tfr(i_c),X(i_c),P,C );
    [SE_ts(i_c).mean,SE_ts(i_c).sd] = get_trial_means(SE_tfr(i_c).trials.TrialSummary,C.t_win.n_sec);
    SE_ts(i_c) = get_mean_ts_feats(SE_ts(i_c),P.feats);
end

%% Run analysis and generate standard figures
if P.vis, spectralevents_vis(SE_tfr,X,TFRs,tVec,P.fVec); end

end

%% helper functions
function [means,sds] = get_trial_means(TS,n_sec)
feats{1} = {'meaneventpower','eventnumber','meaneventduration',     ...
            'meaneventFspan','meanpower'};
feats{2} = {'power','rate','duration','Fspan','trial_power'}; 
n_f = length(feats{1});
class_list = [TS.classLabels]; classes = unique(class_list); 
for i_f = 1:n_f
    feat = feats{2}{i_f}; all_feats = [TS.(feats{1}{i_f})];
    if strcmp(feats{1}{i_f},'eventnumber')
        all_feats = all_feats/n_sec;
    end
    for i_c = classes
        class_feats = all_feats(class_list==i_c);
        means(i_c).(feat) = mean(class_feats,2,'omitnan');
        sds(i_c).(feat)   = std(class_feats,0,2,'omitnan');
    end
end
end

function SE = get_mean_ts_feats(SE,feats)
n_feat = length(feats);
for i_f = 1:n_feat
    try
    feat = feats{i_f}; flist = cell2mat({SE.events.(feat)}');
    catch
        disp(':(')
    end
    SE.mean.(feat) = mean(flist,1,'omitnan');
    SE.sd.(feat) = std(flist,0,1,'omitnan');
end
end