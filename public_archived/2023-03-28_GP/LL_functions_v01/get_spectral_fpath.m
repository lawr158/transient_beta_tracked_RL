% Author: Gerardo Parra, Last updated: 2023-03-28
function fpath = get_spectral_fpath(P,C,chan)
    fname = [int2str(P.method) '_' get_chan_str({},chan) '_' C.t_win.str];
    fname = [fname P.tag '.mat'];
    fpath = fullfile(P.path.out,fname);
end