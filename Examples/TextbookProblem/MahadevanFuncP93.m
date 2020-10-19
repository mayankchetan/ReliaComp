function val = MahadevanFuncP93(x)
% Definition of the performance function for example on P93[]
% ***********************************************************************
% *    Created by Mayank Chetan as a Part of MECH6338 Course at UTD     *
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