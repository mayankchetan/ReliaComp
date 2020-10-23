function [val]=BladeModelTANA(x,varargin)
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
%   This function is used calculate the value of the performance function
%   based on the TANA approximation developed for the MECH6338 project.
%
%   Input:
%           x         = the value about which to evaluate the the function
%                           [(n,1) vector]
%           varargin  = alternative and faster method for calculating the
%                       functional approximation by not reading the MAT 
%                       file in each call. Reduces time for MCS.
%   Output:
%           val       = The value of the performance function at x
%                           [Scalar]
%
%**************************************************************************
if nargin <2
    load('Tana.mat')
else
    TanaAprx = varargin{1};
    secondPoint = varargin{2};
end

% Max deflection allowed
MaxDef = 0.8;

% Calculating the deflection
rR = TanaAprx{3};
Deflection = TanaAprx{1} + (1/rR) .* sum ( (secondPoint.^(1-rR)) .* (x.^rR - secondPoint.^rR) .* TanaAprx{2} );

% Performce function value
val = MaxDef - Deflection;

end