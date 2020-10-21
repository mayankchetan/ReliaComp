function rel = findReliabilityMCS(app)
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

model = app.DataStruct.model;

% Setting up the Covariance matrix.
for idx = 1:length(model.varsMeans)
    for jdx = 1:length(model.varsMeans)
         model.C(idx,jdx) = model.Cprime(idx,jdx) .* model.varSD(idx) .* model.varSD(jdx);
    end
end
% Calculating the Lower triangle matrix [L] using Cholesky decomposition.
model.L = chol(model.C,'lower');

%% Calling the reliability using MCS

for idx = 1:model.MCSamples
    % Sampling the correlated non-normal random variables [1]
    % {X} = [L] * {y} + Mean(X)
    sampleX = model.L * randn (length(model.varsMeans),1) + model.varsMeans;
    g(idx) = feval(model.funcName,sampleX);
    

    
    
    if mod(idx,10000) == 0
        OutputString = sprintf("Finished %d samples of %d \n",idx,model.MCSamples);
        app.UpdateOutputWindow(OutputString,1)
        if CheckKillSwitch(app)
            OutputString = sprintf("Killing Simulation based on user request! \n");
            app.UpdateOutputWindow(OutputString,1)
            break;
        end
    end
    
end

%Finding the number of failures
Nf = sum (g<=0);
rel.Pf = Nf/model.MCSamples;
rel.R = 1 - rel.Pf;

%Calculation of COV and 95% confidence interval
COV = sqrt(rel.R*rel.Pf/model.MCSamples)/rel.Pf;
rel.ConfUp = 1 - (rel.Pf - 1.96 .* (COV.*rel.Pf));
rel.ConfDown = 1 - (rel.Pf + 1.96 .* (COV.*rel.Pf));

end