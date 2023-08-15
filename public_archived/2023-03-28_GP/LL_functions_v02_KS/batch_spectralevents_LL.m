% Author: Gerardo Parra, Last Updated: 2023-03-28
function P = batch_spectralevents_LL(P,C)
%% BATCH_SPECTRALEVENTS_LL: runs spectral events on converted BEAPP batch
% Inputs
% Outputs



%% load batch data and prepare outputs
% get directory converted BEAPP data
P.file.dir = dir( fullfile(C.path.out,'*.mat') ); 
P.file.n = length(P.file.dir);
% set up output directory
P.path.out = get_spectral_fpath(P,C); 
if ~exist(P.path.out,'dir'), mkdir(P.path.out); end

%% iterate through files and run toolbox
for i_f = 1:P.file.n
    F = P.file.dir(i_f);
    % load file data
    load(fullfile(F.folder,F.name),'X')

    % TODO: FOM analysis

    % run spectral events toolbox on data
    [SE_ts,~,TFR] = spectralevents_LL(P,C,X); 

    % TODO: timeseries event detection

    % save toolbox outputs
    save(fullfile(P.path.out,F.name),'SE_ts','X','P','C','TFR')
end

end