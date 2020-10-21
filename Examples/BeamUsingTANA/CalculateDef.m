function TotalDisp = CalculateDef(X)
% Function to calculate the Deflection to blade tip
% ***********************************************************************
% *    Created by Mayank Chetan as a Part of MECH6338 Course at UTD     *
% ***********************************************************************
%   This function is used to calculate the deflection of the loaded wind
%   turbine blade using the provided directional moduli. This function
%   uses a Euler-Bernoulli Beam finite element formulation to calculate the
%   deflection. 
%   
%   A modified version of txt2mat[1] is used to efficiently read large 
%   ASCII data from text files.
%   
%   Input:
%           X         = the value about which to evaluate the the function
%                           [(n,1) vector]
%
%   Output:
%           TotalDisp = The deflection of the blade at the trailing edge
%                       tip. [Scalar]
%
%
% [1] Andres (2019). txt2mat 
%     (https://www.mathworks.com/matlabcentral/fileexchange/18430-txt2mat), 
%     MATLAB Central File Exchange. Retrieved May 3, 2019.
%**************************************************************************
%% Assigning the variables to a more intuitive form.
E11 = X(1);
E22 = X(2);
G12 = X(3);

E33 = 14330e6 .* 0.1;

%% Gathering data from the blades distributed properties. (cross-sectional inertia values)
DistProp=txt2mat('Blade_Distributed_Properties.txt',3); % [1]

%% Processing the simplified beam model and running the static analysis

nele=length(DistProp)-1;

bladeLength=20.8692;

for element=1:nele
    
    % Element Length (m)
    prop(element,1) = (DistProp(element+1,1)-DistProp(element,1))*bladeLength;
    
    % Element Mass Density (kg/m) or rho*A
    prop(element,2) = mean([DistProp(element+1,4) DistProp(element+1,4)]);
    
    % Element Flap-wise Stiffness EI
    prop(element,3) = mean([DistProp(element+1,5) DistProp(element+1,5)]).*E11;
    
    % Element Edgewise Stiffness EI
    prop(element,4) = mean([DistProp(element+1,6) DistProp(element+1,6)]).*E22;
    
    % Element Torsional Stiffness EI
    prop(element,5) = mean([DistProp(element+1,7) DistProp(element+1,7)]).*G12;
    
    % Element EA
    prop(element,6) = mean([DistProp(element+1,8) DistProp(element+1,8)]).*E33;
    
    % Element Twist
    prop(element,7) = mean([DistProp(element+1,3) DistProp(element+1,3)]);
    
end

% Fetching the Global Stiffness Matrix
[~,K] = ortho_beam_fem(nele,prop);

% Applying the cantilever BC
K(1:6,:) = 0;
K(:,1:6) = 0;
K(1:6,1:6) = eye(6);

force = zeros((nele+1)*6,1);

% Load on nodes as defined by the NREL ground testing for t = 200s
loadNodes = [23 26 30];
loadVals = [8.681465476000000e+03,2.858377809000000e+03,7.093107165000000e+03]; %@ t=200s on Max Load condition

% appltying in the z or Flap-wise direction
force(loadNodes.*6-3) = loadVals;

% appltying in the y or Flap-wise direction
force(loadNodes.*6-4) = loadVals;

% applying in the Twist direction to account for elastic acis offset
force(loadNodes.*6-2) = loadVals(:) .* DistProp(loadNodes,17).*5;

% Calculating the Displacements
Disp=K\force;

% Arranging the displacements in a usable form
displacementSimple=[Disp(1:6:end) Disp(2:6:end) Disp(3:6:end)];
displacementSimpleRot=[Disp(4:6:end) Disp(5:6:end) Disp(6:6:end)];

% Disp of last node and trailing edge!!
LastChord = 0.360086064; % of the tip

ZDispRot = LastChord .* sin(displacementSimpleRot(end,1));
yDispRot = LastChord .* (cos(displacementSimpleRot(end,1))-1);

TotalDisp = sqrt((displacementSimple(end,3)+ZDispRot).^2+(displacementSimple(end,2)+yDispRot).^2);


% Code to visualize the displacement, set below if to 1
if 0
    f=figure;
    plot(DistProp(:,1)*bladeLength, displacementSimple(:,3)*1,'g','LineWidth',1.5);
    hold on
    plot(DistProp(:,1)*bladeLength, displacementSimpleRot(:,[1:3])*1,'r','LineWidth',1.5);
    axis([0 24 0 1])
    xlabel('blade length');
    ylabel('flap-wise displacement')
    grid on
    grid minor
end
end

%% All Internal functions
function [uNorm] = calculateNorm(uPrev,u)
% This function calculates the relative norm between two vectors
%
%    input:
%    uPrev   = displacement vector from previous step
%    u       = displacement vector from current step
%
%    output:
%    uNomrm  = calculated relative norm of u - uPrev

len = length(u);
numerator = 0;
denom = 0;
for i=1:len
    numerator = numerator + (u(i)-uPrev(i))^2;
    denom = denom + u(i)^2;
end

uNorm = sqrt(numerator/denom);
end