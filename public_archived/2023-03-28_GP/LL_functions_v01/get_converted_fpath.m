% Author: Gerardo Parra, Last updated: 2023-03-28
function fpath = get_converted_fpath(C,chan)
    fname = [get_chan_str({},chan) '_' C.t_win.str '.mat'];
    fpath = fullfile(C.path.out,'inputs',fname);
end