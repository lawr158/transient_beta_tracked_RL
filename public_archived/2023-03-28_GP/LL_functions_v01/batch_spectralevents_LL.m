% Author: Gerardo Parra, Last Updated: 2023-03-28
function P = batch_spectralevents_LL(P,C)
%% BATCH_SPECTRAL: runs spectral events on converted BEAPP data
% Inputs
% Outputs

%%
fRange = ['_' get_range_str(P.band.range,'Hz')];
P.path.out = fullfile(C.path.out,[P.band.name fRange]);
if ~exist(P.path.out,'dir'), mkdir(P.path.out); end

for i_c = 1:C.chan.n
    chan = C.chan.list(i_c);
    % load channel data
    load(get_converted_fpath(C,chan),'X','T')

    % TODO: FOM analysis

    % run spectral events toolbox on data
    [SE,TFR,~] = spectralevents_LL( P.band.range,P.fVec,C.Fs,P.method,  ...
                                    P.vis,X,T,P.ncyc ); 

    % TODO: timeseries event detection

    % save toolbox outputs
    save(get_spectral_fpath(P,C,chan),'SE','TFR','P')
end

end