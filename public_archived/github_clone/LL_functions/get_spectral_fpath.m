% Author: Gerardo Parra, Last updated: 2023-03-29
function fpath = get_spectral_fpath(P,C)
%% GET_SPECTRAL_FPATH: get filepath for spectral events toolbox
% Inputs
%   - P - struct with spectral events batch parameters
%   - C - struct with BEAPP conversion parameters
% Outputs
%   - fpath - string with folder path
%%
% analysis frequency band
band = [P.band.name '_' get_range_str(P.band.range,'Hz')];
% time window and batch tag
t_win = [C.t_win.str P.tag];
% event-finding method
method = ['method' int2str(P.method)];
% concatenate information in folder name
fname = [band '_' t_win '_' method];
% create path from batch path and folder name
fpath = fullfile(C.path.in,['spectral_events' C.tag],fname);

end