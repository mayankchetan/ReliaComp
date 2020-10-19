% ***********************************************************************
% *    Created by Mayank Chetan as a Part of MECH6338 Course at UTD     *
% ***********************************************************************
% This runs a Monte Carlo Simulations (MCS) approximation to find the 
% reliability related to the MECH6338 project. The function uses the 
% Cholesky decomposition method[1]. The functional approximation MAT file 
% obtained by running the "FitTANA.m" script is used here.
% 
%
% [1] Haldar, Achintya, and Sankaran Mahadevan. Reliability assessment 
%     using stochastic finite element analysis. John Wiley & Sons, 2000.
%**************************************************************************
%% Setup for reliability analysis
tic
% adds path to all functions used by this script.
addpath('utils');

% the input is defined in terms of a structure called 'model'
model.funcName = 'PerfFunc'; % Name of the performance function
%                     [E11 , E22 , G12]
model.varsMeans = [35910e6;14330e6;3510e6];  % Means for the r.v's
model.varsCOV = [0.05; 0.05; 0.05]; % COV for the r.v's
model.gradFlag = 1; % Flag = 1 for numerical gradients (Central Difference)
model.gradDelta = 1e-6; % numerical gradient perturbation in %
model.BetaDiff = 1e-9; % Convergence criteria for Beta
model.AlphaDiff = 1e-9; % Convergence criteria for Alpha, >1 to turn off

model.varSD = model.varsMeans.*model.varsCOV; % Calc of SD's

% Distribution Type 0 = Normal, 1 = Log Normal, 2 = Uniform
model.Types = [1;1;1];

% Number of samples for MCS run
SamplesN = 1e6;

% Assuming no corrections to the Correlation Coefficient
rho_E11_E22 = 0.1;
rho_E11_G12 = 0.15;
rho_E22_G12 = 0.05;

model.Cprime = [1               rho_E11_E22     rho_E11_G12;
                rho_E11_E22     1               rho_E22_G12;
                rho_E11_G12     rho_E22_G12     1          ];

% Setting up the Covariance matrix.
COV_E11_E22 = rho_E11_E22 .* model.varSD(1) .* model.varSD(2);
COV_E11_G12 = rho_E11_G12 .* model.varSD(1) .* model.varSD(3);
COV_E22_G12 = rho_E22_G12 .* model.varSD(2) .* model.varSD(3);

model.C  =[ model.varSD(1)^2    COV_E11_E22         COV_E11_G12
            COV_E11_E22         model.varSD(2)^2    COV_E22_G12
            COV_E11_G12         COV_E22_G12         model.varSD(3)^2];
        
% Calculating the Lower triangle matrix [L] using Cholesky decomposition.
model.L = chol(model.C,'lower');

%% Calling the reliability using MCS
load('Tana.mat')

for ii = 1:SamplesN
    
    % Sampling the correlated non-normal random variables [1]
    % {X} = [L] * {y} + Mean(X)
    sampleX = model.L * randn (3,1) + model.varsMeans;
    g(ii) = PerfFunc(sampleX,TanaAprx,secondPoint);
end

%Finding the number of failures
Nf = sum (g<=0);
Pf = Nf/SamplesN;
R = 1 - Pf;

%Calculation of COV and 95% confidence interval
COV = sqrt(R*Pf/SamplesN)/Pf;
ConfUp = 1 - (Pf - 1.96 .* (COV.*Pf));
ConfDown = 1 - (Pf + 1.96 .* (COV.*Pf));

fprintf('The reliability for %2.2E MCS runs is %2.6f with 95%% CI (%2.6f,%2.6f)  \n',SamplesN,R,ConfDown,ConfUp)

toc
