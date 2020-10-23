function [val]=DeflectionFunc(x)
% Definition of the performance function
% ***********************************************************************
% Copyright (C) 2020  Mayank Chetan
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
% ***********************************************************************
%   This function is used to define the performance function g. It
%   calculates the value of it for the point passed. Currently this act
%   like an anonymous function.
%   
%   Input:
%           x         = the value about which to evaluate the the function
%                           [(n,1) vector]
%
%   Output:
%           val       = The value of the performance function at x
%                           [Scalar]
%
%**************************************************************************


% Value of Norm used for fitting the TANA aproximation
val = CalculateDef(x);

end