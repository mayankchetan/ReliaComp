% ***********************************************************************
% *    Created by Mayank Chetan as a Part of MECH6338 Course at UTD     *
% ***********************************************************************
% This runs a Hasofer-Lind-Rackwitz-Fiessler (HL-RF) algorithm to find the 
% reliability related example on P93[1]. The function uses the 
% Eigen value method[1] to calculate the reliability analytically. 
%
% 
% 
%
% [1] Haldar, Achintya, and Sankaran Mahadevan. Reliability assessment 
%     using stochastic finite element analysis. John Wiley & Sons, 2000.
%**************************************************************************
%% Setup for reliability analysis

% the input is defined in terms of a structure called 'model'
addpath('utils');
model.funcName = 'MahadevanFuncP93'; % Name of the performance function
% [E11 , E22 , G12]
model.varsMeans = [38;54];  % Means for the r.v's
model.varsCOV = [0.1;0.05]; % COV for the r.v's
model.gradFlag = 1; % Flag = 0 for analytical, 1 for numerical gradients
model.gradDelta = 1e-6; % numerical gradient perturbation in %
model.BetaDiff = 1e-9; % Convergence criteria for Beta
model.AlphaDiff = 1e-6; % Convergence criteria for Beta

model.varSD = model.varsMeans.*model.varsCOV; % Calc of SD's

% Dist Type 0 = Normal, 1 = Log Normal
model.Types = [1;0];

% Assuming no corrections to the Correlation Coefficient

rho_E11_E22 = 0.1;
rho_E11_G12 = 0.15;
rho_E22_G12 = 0.05;

model.Cprime = [1       0.3
                0.3     1];
            
model.LNSD = sqrt(log(1 + model.varsCOV.^2));
model.LNmean = log(model.varsMeans) - 0.5 .* model.LNSD .^2;
%% Calling the reliability calculator

% Returned value is a structure containing data for all iterations.
reliability = findReliability(model);

fprintf('The final Beta = %2.8f, with a reliability of %2.8f\n',reliability(end).beta,reliability(end).R)

