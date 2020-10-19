function rel = HLRFAlgo(model,rel,app)
% function to calculate the reliability using HLRF method
% ***********************************************************************
% *    Created by Mayank Chetan as a Part of MECH6338 Course at UTD     *
% ***********************************************************************
%   This function calculates the reliability of the provided function using
%   the HL-RF method[1]. most calculations are vectorized
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
%           rel        =  the reliability structure passed tp this
%                           function. [This is defined below] [ix1 struct]
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
%           rel.ADelta =  change in alpha within each iteration of beta.
%                           [Scalar]
%           rel.U =  The value of U's for the ith iteration.
%                           [(n,1) vector]
%
% [1] Haldar, Achintya, and Sankaran Mahadevan. Reliability assessment 
%     using stochastic finite element analysis. John Wiley & Sons, 2000.
%**************************************************************************

% Finding the current iteration number
k = length(rel)+1;

%% Alpha Convergence for more stability, limited to 100 iterations.
alphaCount = 1;
while (rel(k-1).ADelta >= model.AlphaDiff) && alphaCount <=100
    
    OutputString = sprintf("Beta Iter = %d, Alpha Iter = %d \n",k,alphaCount);
    app.UpdateOutputWindow(OutputString,1)
    
    %fprintf('Beta Iter = %d, Alpha Iter = %d \n',k,alphaCount)
    
    % Finding the Equivalent Normal Mean and SD depending on Distribution.
    % This runs through the different random variables.
    for ii=1:length(rel(k-1).xStar)
        
        % Log-normal Distribution
        if model.Types(ii)==1
            DistType = 'Lognormal';
            
            xStar = rel(k-1).xStar(ii);
            
            LNMean  = model.LNmean(ii);
            LNsd    = model.LNSD(ii);
            
            LNCDFVal = cdf(DistType,xStar,LNMean,LNsd);
            LNPDFVal = pdf(DistType,xStar,LNMean,LNsd);
            
            equNormSD(ii) = pdf('Normal',norminv(LNCDFVal),0,1)/LNPDFVal;
            
            equNormMean(ii) = xStar - norminv(LNCDFVal) .* equNormSD(ii);
            
        % Uniform Distribution
        elseif model.Types(ii)==2
            DistType = 'Uniform';
            
            xStar = rel(k-1).xStar(ii);
            
            Upper  = model.varsMeans(ii) + model.varSD(ii);
            Lower    = model.varsMeans(ii) - model.varSD(ii);
            
            % A saturation function is added to the uniform distribution to
            % avoid the value of Xstar to be beyond the bounds resulting in
            % a p(x)->0 hence causing EquNormSD -> Inf
            
            % Saturating the value of xStar to 1% within bounds
            xStar = max (Lower*1.01, min(xStar,Upper*0.99));
            
            UCDFVal = cdf(DistType,xStar,Lower,Upper);
            UPDFVal = pdf(DistType,xStar,Lower,Upper);
            
            equNormSD(ii) = pdf('Normal',norminv(UCDFVal),0,1)/UPDFVal;
            
            equNormMean(ii) = xStar - norminv(UCDFVal) .* equNormSD(ii);
            
        % Normal Distribution
        else
            DistType = 'Normal';
            
            xStar = rel(k-1).xStar(ii);
            equNormSD(ii) = model.varSD(ii);
            equNormMean(ii) = model.varsMeans(ii);
            
        end
    end
    
    % Diagonalizing the Equ Normal SD values
    EquNormSD_Diag = diag(equNormSD);
    
    % Finding the Eigen Val and Normalized Eigen Vectors for the Cprime Matrix
    [V,D] = eig(model.Cprime);
    
    % Finding the Reduced Form of the original Variables
    xStar = rel(k-1).xStar;
    xStar_Reduced = (xStar(:) - equNormMean(:)) ./ equNormSD(:); % Making it a vector
    
    % Calculating the Uncorrelated Normal Variables
    y_Vec = V' * xStar_Reduced;
    
    % Symbolically obtaining the value and gradients
    A = sym('y%d',[1 ii]);
    
    SymsXDomain = EquNormSD_Diag * V * A' + equNormMean(:);
    g_y = feval(model.funcName,SymsXDomain);
    dg_gy = gradient(g_y,A);
    
    dg_gy_Val = double(subs(dg_gy, A', y_Vec));
    
    
    % getting the SD for g using the equivalent normal variables.
    SDG = sqrt((dg_gy_Val.^2)'*diag(D));
    
    % Saving Previous Alpha
    Temp = rel(k-1).alpha;
    
    % Calculating Alpha's ***POSITIVE ALPHA HERE!!! CAREFUL!!!!***
    rel(k-1).alpha =(dg_gy_Val.*sqrt(diag(D)))./SDG';
    
    % Change in alpha, used for alpha convergence
    rel(k-1).ADelta =max(abs(Temp - rel(k-1).alpha));
    
    % Ystar = - Alpha * [Sigma] * Beta
    Ystar = - rel(k-1).alpha .* sqrt(diag(D)) .* rel(k-1).beta;
    
    % New Xstar values in Non-Normal, Correlated domain.
    rel(k-1).xStar = EquNormSD_Diag * V * Ystar + equNormMean(:);
    
    alphaCount = alphaCount+1;
end

% ALPHA FINALLY CONVERGED!!!!!!
rel(k).xStar = rel(k-1).xStar;
rel(k).alpha = rel(k-1).alpha;

% Calculating a Beta for which g(Xstar) = 0
syms BetaNew

YstarNew_Sym = - rel(k-1).alpha .* sqrt(diag(D)) .* BetaNew;
XstarNew_Sym = EquNormSD_Diag * V * YstarNew_Sym + equNormMean(:);
g_x_Sym = feval(model.funcName,XstarNew_Sym);

tempFuncBeta = matlabFunction(g_x_Sym);
options = optimoptions('fsolve','Display','off','FunctionTolerance',1e-12);
[rel(k).beta,error] = fsolve(tempFuncBeta,rel(k-1).beta*0.1,options);

% Calculating Reliability
rel(k).R = normcdf(rel(k).beta,0,1,[]);
rel(k).Pf = 1-rel(k).R;

% Calculating the change in beta wrt previous iteration.
rel(k).BDelta = abs(rel(k).beta-rel(k-1).beta);

% Setting change in alpha for next iteration as 1.
rel(k).ADelta = 1;
end