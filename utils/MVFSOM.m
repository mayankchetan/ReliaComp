function rel = MVFSOM(model)
% function to calculate the reliability using MVFSOM Method
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
%   This function calculates the reliability of the provided function using
%   the MVFOSM method
%   
%   Input:
%           model      =  the problem structure passed to this function.
%                           [This is defined below]
%           model.funcName  = name of the performance function passed to 
%                               this solver [String]
%           model.varsMeans = the Mean values for the n random variables 
%                               [(n,1) vector]
%           model.varsSD    = the standard deviation values for the n 
%                               random variables [(n,1) vector]
%           model.gradFlag  = Flag for analytical (0) or numerical (1)
%                               gradient calculation method [flag]
%           model.gradDelta = numerical gradient perturbation in percentage
%                               of the r.v [scalar]
%           model.BetaDiff  = Convergence criteria for 2 consecutive beta
%                               calculations [scalar]
%
%   Output:
%           rel        =  the reliability structure passed from this 
%                           function. [This is defined below] [ix1 struct]
%           rel.xStar  =  The value of Xstar for the ith iteration. 
%                           [(n,1) vector]
%           rel.grad   =  The value of gradient at the Xstar point. 
%                           [(n,1) vector]
%           rel.GVal  =  The value of the performance function at Xstar
%                           [Scalar]
%           rel.SDG   =  The value of 'SD' of performance function.
%                           sqrt(sum(grad .* SD(rv's))) [Scalar]
%           rel.beta  =  The value of beta for the ith iteration. 
%                           [Scalar]
%           rel.R     =  The value of reliability for the ith iteration. 
%                           [Scalar]
%           rel.Pf    =  The value of probability of failure for the ith 
%                           iteration.[(n,1) vector]
%           rel.alpha =  The value of alpha's at the previous Xstar point. 
%                           [(n,1) vector]
%           rel.BDelta=  change in beta from the previous iteration. 
%                           [Scalar]
%           rel.alpha =  The value of U's for the ith iteration. 
%                           [(n,1) vector]
%
%**************************************************************************
%% Setting up the variables
rel(1).xStar = model.varsMeans;
SDx = model.varSD;

%% Calculating the gradients
[rel(1).grad,rel(1).GVal]=FindGrads(model.funcName,...
                          rel(1).xStar,model.gradDelta);

%% MVFSOM Calculation

%finding the SD
rel(1).SDG = sqrt( sum( (rel(1).grad.*SDx).^2 ) );

% Calculation of beta, R, Pf
rel(1).beta = rel(1).GVal/rel(1).SDG;
rel(1).R = normcdf(rel(1).beta,0,1,[]);
rel(1).Pf = 1-rel(1).R;

% Setting alpha, BDelta, U to zero.
rel(1).alpha = zeros(length(rel(1).xStar),1);
rel(1).BDelta = 0;
rel(1).ADelta = 1;
rel(1).U = zeros(length(rel(1).xStar),1);
end