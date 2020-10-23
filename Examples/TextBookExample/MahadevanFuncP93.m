function val = MahadevanFuncP93(x)
% Definition of the performance function for example on P93[1]
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
%   This function is used calculate the value of the performance function
%   based on the example on page 93 of [1].
%
%   Input:
%           x         = the value about which to evaluate the the function
%                           [(n,1) vector]
%   Output:
%           val       = The value of the performance function at x
%                           [Scalar]
%
% [1] Haldar, Achintya, and Sankaran Mahadevan. Reliability assessment 
%     using stochastic finite element analysis. John Wiley & Sons, 2000.
%**************************************************************************

Fy = x(1);
Z = x(2);

val = Fy .* Z - 1140;

end