function time_str = get_time_win_str(time_win)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
time_str = strrep(int2str(time_win),'   ',' ');
time_str = strrep(time_str,'  ',' ');
time_str = strrep(time_str,' ','to');
time_str = [time_str 'ms'];
end

