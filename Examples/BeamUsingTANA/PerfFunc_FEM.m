function [val]=PerfFunc_FEM(x)
% Definition of the performance function
% ***********************************************************************
% *    Created by Mayank Chetan as a Part of MECH6338 Course at UTD     *
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
% Max deflection allowed
MaxDef = 0.8;

% Calculating the deflection
Deflection = CalculateDef(x);

% Performce function value
val = MaxDef - Deflection;

end