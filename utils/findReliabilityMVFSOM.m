function rel = findReliabilityMVFSOM(app)
% function to calculate the reliability using MVFOSM
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
% ************************************************************************
%   This function calculates the reliability of the provided function using
%   the MVFOSM methods for the first iteration assuming that the variables 
%   are independent and normally distributed to obtain a Beta approx.
%   The Hasofer-Lind-Rackwitz-Fiessler (HL-RF) method [1] is used for the 
%   subsequent iterations (limited to 10,000 iterations).The underlying 
%   functions are not limited to any specific performance function. the
%   user can define a custom function and pass it to this solver.
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
%           model.AlphaDiff = Convergence criteria for 2 consecutive alpha
%                               calculations [scalar]
%           model.LNMeans   = the LN Mean values for the n random variables
%                               [(n,1) vector]
%           model.LNSD      = the LN standard deviation values for the n
%                               random variables [(n,1) vector]
%           model.CPrime    = the correlation coefficient for the random 
%                               variables [(n,n) matrix]
%           model.types     = the distribution types for the random 
%                               variables, Normal, LN or Uniform [(n,1)Vec]
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
%           rel.ADelta =  change in alpha within each iteration of beata.
%                           [Scalar]
%           rel.U =  The value of U's for the ith iteration.
%                           [(n,1) vector]
%
% [1] Haldar, Achintya, and Sankaran Mahadevan. Reliability assessment 
%     using stochastic finite element analysis. John Wiley & Sons, 2000.
%**************************************************************************
idx=1;

model = app.DataStruct.model;

% finding reliability based on MVFOSM assuming Uncorrelated Normal Dist
rel(idx) = MVFSOM(model);
       
end