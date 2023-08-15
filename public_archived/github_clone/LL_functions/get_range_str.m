% Author: Gerardo Parra, Last updated: 2023-03-29
function range_str = get_range_str(range,varargin)
%%  GET_RANGE_STR: returns string of numeric range
% Inputs:
%   - range - 1x2 numeric array with range to get string for
%   - unit  - optional, char with range units to append at end of string
% Outputs
%   - range_str - string with the format [range(1) 'to' range(2) units]
%%
% set unit
if nargin > 1, unit = varargin{1}; else, unit = ''; end
range_str = [sprintf('%dto%d',range(1),range(2)) unit];

end

