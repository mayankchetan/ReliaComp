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
% This fits a TANA approximation for the performance function related to the
% MECH6338 project. The constants related to the approximation is stored in
% a .mat file that is read in by the other codes developed for calculating
% the reliability.
%**************************************************************************
close all
clear

% adds path to all functions used by this script.
addpath('../BladeModelFEM');
addpath('../../utils');

% Expansion Point defined
E11 = 35910e6;  %Pa
E22 = 14330e6;  %Pa
G12 = 3510e6;   %Pa

expansionPoint = [E11 E22 G12]';

% scaling variable to validate the approximation
t = 0.6:0.05:1.2;

%% Calculating values for the approximations

% Calculates the TANA approximation at the expansion point
[TanaAprx] = findFuncApprox(expansionPoint,4);

for ii = 1:length(t)
    
    % Running through the different points to evaluate the accuracy of
    % approximation.
    evalPoint = expansionPoint .* t(ii);
    
    % TANA Approx points
    rR = TanaAprx{3};
    secondPoint = expansionPoint .* 0.8;
    TanaFuncVals(ii) = TanaAprx{1} + (1/rR) .* sum ( (secondPoint.^(1-rR)) .* (evalPoint.^rR - secondPoint.^rR) .* TanaAprx{2} );

    ExactVals(ii) = DeflectionFunc(evalPoint);
    
end

figure(455)

plot(t,ExactVals,'-pr');
hold on
plot(t,TanaFuncVals,'-*k');
ylabel('Tip Displacement (m)')

yyaxis right
plot(t,(ExactVals-TanaFuncVals)./ExactVals*100,'--k');

legend('Exact','TANA','% Error for TANA','Location','NorthWest')
grid on
grid minor
xlabel('Magnitude (t) of Direction vector (-)')
ylabel('% Error for approximation (%)')

title('Displacement of tip (m)')

% Saving the TANA approximation as a MAT file to be used by other routines
save('Tana','TanaAprx','secondPoint')

copyfile('Tana.mat','../BladeModelTANA');