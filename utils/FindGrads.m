function [grad,val]=FindGrads(funcName,x,delta)
% function to calculate the gradients for the function passed at point
% ***********************************************************************
% *    Created by Mayank Chetan as a Part of MECH6338 Course at UTD     *
% ***********************************************************************
%   This function calculates the gradients for the passed function at the
%   defined point using analytical or central difference method.
%
%   Input:
%           funcName  = name of the performance function passed to
%                           this solver [String]
%           x         = the value about which to evaluate the gradients
%                           [(n,1) vector]
%           delta     = numerical gradient perturbation in percentage
%                           of the r.v [scalar]
%
%   Output:
%           grad      = The value of gradients at the point x.
%                           [(n,1) vector]
%           val       = The value of the performance function at Xstar
%                           [Scalar]
%
%**************************************************************************

% Finding the value of the function at the point x
[val]=feval(funcName,x);

% Central Difference Gradients
xGrad = x;
xUp = x;
xDown = x;

% Moving through the different r.v's
for ii = 1:length(x)
    
    if xGrad(ii) == 0
        xGrad(ii) = 1e6;
    end
    
    DeltaX=xGrad(ii)*delta;
    xUp(ii) = xGrad(ii) + DeltaX/2;
    xDown(ii) = xGrad(ii) - DeltaX/2;
    valUp=feval(funcName,xUp);
    valDown=feval(funcName,xDown);
    
    grad(ii) = (valUp-valDown)/(DeltaX);
    xUp = x;
    xDown = x;
end

% Converting into a column vector
grad=grad(:);

end