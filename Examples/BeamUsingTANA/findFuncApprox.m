function varargout=findFuncApprox(expansionPoint,approxMethod)
% Function to calculate the constants for functional approximation
% ***********************************************************************
% *    Created by Mayank Chetan as a Part of MECH6338 Course at UTD     *
% ***********************************************************************
%   This function is used calculate the constants used in types of
%   functional approximations methods described below.
%
%   Approximation methods:
%   4. TANA
%
%   Input:
%           expansionPoint  = Vector containing the expansion point for the
%                           approximation [(n,1) vector]
%           approxMethod    = Method used for the Approximation
%                           [scalar]
%   Output:
%           varargout = Constants related to the selected approximation
%                       method, two cells represent the radial and
%                       tangential stress respectively [1x2 cell]
%
%**************************************************************************

%% Assigning the Values for the approximations

switch(approxMethod)
    
    case 4 % TANA Approx
        
        % Value of the function at the expansion point
        GRX1 = DeflectionFunc(expansionPoint);
        
        % 20% lower from expansion point
        secondPoint = expansionPoint .*0.8; 

        % Function value g(X2) & Gradients dg/dxi(X2)
        [varargout{1}{2},varargout{1}{1}] = FindGrads('DeflectionFunc',secondPoint,0.01);
        
        % Finding the value of R
        % Setting up equation
        syms R
        
        TANAr = (varargout{1}{1} + (1/R) .* sum( (secondPoint.^(1-R)).* (expansionPoint.^R - secondPoint.^R) .* varargout{1}{2})) - GRX1;
        tempFuncr = matlabFunction(TANAr);
        options = optimoptions('fsolve','Display','off','FunctionTolerance',1e-9);
        [varargout{1}{3},fval] = fsolve(tempFuncr,1);
        
        if abs(fval) > 1e-5
            warning('Error for TANA approximation is high -> %d',fval)
        end
        
        
    case 5 % TANA1 Approx
        
        % Function value g(X1) & Gradients dg/dxi(X1)
        [varargout{1}{2},varargout{1}{1}] = FindGrads('DeflectionFunc',expansionPoint,0.01);
       
        secondPoint = expansionPoint .* 0.7; % 20% lower from expansion point
        
        syms p
        pR = 5*ones(length(expansionPoint),1);
        
        [Grad_SigrSecondPt,GRX2] = FindGrads('DeflectionFunc',secondPoint,0.01); % This will be gradients at the second point
        
        for ii = 1:length(pR)
            
            TANA1r = ((secondPoint(ii)/expansionPoint(ii))^(p-1))*varargout{1}{2}(ii) - Grad_SigrSecondPt(ii);
            tempFuncr = matlabFunction(TANA1r);
            options = optimoptions('fsolve','FunctionTolerance',1e-20);
            [pR(ii),fval] = fsolve(tempFuncr,pR(ii),options);
        end
        
        varargout{1}{3} = pR;
        
        ErrorR = GRX2 - (varargout{1}{1} + sum((varargout{1}{2}).*((expansionPoint.^(1-pR))./pR).*(secondPoint.^pR - expansionPoint.^pR)));
        
        varargout{1}{4} = ErrorR;
        
        
end
end