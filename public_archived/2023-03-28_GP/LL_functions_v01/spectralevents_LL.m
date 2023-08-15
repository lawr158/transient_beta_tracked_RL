function [SE,TFRs,X] = spectralevents_LL( fBand,fVec,Fs,findMethod, ...
                                            vis,varargin )
%% SPECTRALEVENTS_LL Find and analyze transient spectral events in a 
%   time-series dataset. Spectral events are defined as local maxima above 
%   a power threshold of a specified band in the non-averaged 
%   time-frequency responses (TFR). Modified for data from the Levin Lab.
%   Last updated by: Gerardo Parra, 2023-03-29
%
%  [SE,TFRs,X] = SPECTRALEVENTS(fBand,fVec,Fs,findMethod,vis,X,classLabels)
%   or
%  [SE,TFRs,X] = SPECTRALEVENTS(fBand,fVec,Fs,findMethod,vis,X{1}, ...
%                               classLabels{1},X{2},classLabels{2},...)
%   
%   Returns a structure array of spectral event features (specEv_struct),
%   cell array containing the time-frequency responses (TFRs), and cell 
%   array of all time-series trials (X) for each subject/session within the
%   dataset comparing various experimental conditions or outcome states 
%   corresponding to each trial. Note: this function sets the factors of
%   median threshold (thrFOM) at 6.
%
% Inputs:
%   eventBand - range of frequencies ([Fmin_event Fmax_event]; Hz) over 
%       which above-threshold spectral
%       power events are classified.
%   fVec - frequency vector (Hz) over which the time-frequency response 
%       (TFR) is calculated. Note that this set must fall within the range 
%       of resolvable/alias-free frequency values (i.e. Fmin>=1/(trial 
%       duration), Fmax<=(Nyquist freq)).
%   Fs - sampling frequency (Hz).
%   findMethod - integer value specifying which event-finding method to use
%       (1, 2, or 3). Note that the method specifies how much overlap
%       exists between events. Use 1 to replicate the method used in 
%       et al. eLife 2017.
%   vis - logical value that determines whether to run basic feature 
%       analysis and output standard figures.
%   X{a} - m-by-n matrix (of the a^th subject/session cell in cell array X) 
%       representing the time-series trials of the given subject. m is the 
%       number of timepoints and n is the number of trials. Note that m 
%       timepoints must be uniform across all trials and subjects.
%   T{a} - numeric or logical trial classification labels (of the
%       a^th subject/session cell in cell array classLabels); associates 
%       each trial of the given subject/session to an experimental 
%       condition/outcome/state (e.g., hit or miss, detect or non-detect, 
%       attend-to or attend away). If classLabels{a} is entered as a single
%       value (e.g., 0 or 1), all trials in the a^th subject/session are 
%       associated with that label. Alternatively, classLabels{a} can be 
%       entered as a vector of n elements, each corresponding to a trial 
%       within the a^th subject/session.
%   ncyc - optional, integer number of cycles to use for wavelet transform 
%   thrFOM - optional, integer to use as FOM threshold
%
% Outputs:
%   specEv_struct - array of event feature structures, each corresponding
%       with a subject/session, respectively.
%   TFRs - cell array with each cell containing the time-frequency response 
%       (freq-by-time-by-trial) for a given subject/session.
%   X - cell array with each cell containing the time-series trials for a
%       given subject/session.
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

% set default wavelet cycle number if not included in inputs
if nargin >= 8, ncyc = varargin{3}; else, ncyc = 7; end
% set default Factors of Median threshold (see Shin et al. eLife 2017 for 
%   details concerning this value)
if nargin >= 9, thrFOM = varargin{4}; else, thrFOM = 6; end

% Validate number of time-series (X{1}, X{2},...) and trial class label 
%   (classLabels{1}, classLabels{2},...) inputs
if nargin-5>=2
    % struct array containing time-series trials by subject/session
    X = varargin{1};
    % struct array containing trial classification cells by subject/session
    T = varargin{2}; 
else
    error(['Must input at least one time-series dataset and its '       ...
           'corresponding trial classification label/vector.'])
end

% Validate format of trial class labels and reformat if necessary
n_subj = numel(X); % Number of subjects/sessions
for i_s=1:n_subj
    if ~(isnumeric(T(i_s).class) ||                           ...
            islogical(T{i_s}.class))
        error('Trial classification labels must be numeric or logical.')
    end
    cl_s = T(i_s).class;
    if numel(cl_s)==1
        % Reformat if entered as single values instead of vectors
        cl_s = ones(1,size(X{i_s},2))*cl_s; T(i_s).class = cl_s;
    elseif isequal(size(cl_s),[size(X(i_s).data,2),1])
        T(i_s).class = cl_s';
    elseif ~isequal(size(cl_s),[1,size(X(i_s).data,2)])
        error(['Trial classification labels must be formatted as '      ...
               'either a 1-row array or single value.'])
    end
end
clear cl_s

% Validate fVec input
Fn = Fs/2; % Nyquist frequency
dt = 1/Fs; % Sampling time interval
Fmin = 1/(size(X(1).data,1)*dt); % Minimum resolvable frequency
if fVec(1)<Fmin || fVec(end)>Fn || abs(fVec(2)-fVec(1))<Fmin
    error(['Frequency vector includes values outside the resolvable/'   ...
           'alias-free range.'])
end

% Solve for the time-frequency response (TFR) and run spectral event
%   analysis on each trial within each subject/session
TFRs = struct('id',{X.id},'tfr',[]); % Struct for TFRs across subjects
SE   = struct('id',{X.id},'events',[]);
for i_s=1:n_subj
    TFR = []; % Matrix for storing freq-by-time-by-trial
    for i_t=1:size(X(i_s).data,2)
        % Transformation calculated using a Morlet wavelet (width=7), see 
        %   4DToolbox for function details)
        [TFR_t,tVec,~] = spectralevents_ts2tfr( X(i_s).data(:,i_t),fVec,...
                                                Fs,ncyc ); 
        % Concatenate each trial along the 3rd dimension
        TFR = cat(3,TFR,TFR_t); 
    end
    TFRs(i_s).tfr = TFR; % Append TFR for the given subject
    % Find spectral events
    [SE(i_s).trials,SE(i_s).events,SE(i_s).IEI] =                       ...
        spectralevents_LL_find( findMethod,fBand,thrFOM,tVec,fVec,TFR,  ...
                                T(i_s).class );
end

% Run analysis and generate standard figures
if vis, spectralevents_vis(SE,X,TFRs,tVec,fVec); end

end