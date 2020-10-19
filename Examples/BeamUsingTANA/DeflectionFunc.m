function [val]=DeflectionFunc(x)
% Definition of the performance function
% ***********************************************************************
% *    Created by Mayank Chetan as a Part of MECH6338 Course at UTD     *
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