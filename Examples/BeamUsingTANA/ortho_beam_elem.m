
function [Me,Ke] = ortho_beam_elem(L,MassDen,FlapEI,EdgeEI,GJ,EA,Twist)

% 2 node beam element
% 6 DOF per node
% Based on definition in
% Przemieniecki, J. S. “Theory of Matrix Structural Analysis” . New York: Dover, 1985. Print.
%
% "x" is the beam axis
%
% Me and Ke are element mass and stiffness matrices

Me = zeros(12,12);
% Me(1,1) = 140;  Me(1,7) = 70;
% Me(2,2) = 156;  Me(2,6) = 22*L; Me(2,8) = 54;   Me(2,12) = -13*L;
% Me(3,3) = 156;  Me(3,5) = -22*L;    Me(3,9) = 54;  Me(3,11) = 13*L;
% Me(4,4) = 140*Ip/A; Me(4,10) = 70*Ip/A;
% Me(5,3) = Me(3,5);  Me(5,5) = 4*L^2;    Me(5,9) = -13*L;    Me(5,11) = -3*L^2;
% Me(6,2) = Me(2,6);  Me(6,6) = 4*L^2;    Me(6,8) = 13*L;     Me(6,12) = -3*L^2;
% Me(7,1) = Me(1,7);  Me(7,7) = 140;
% Me(8,2) = Me(2,8);  Me(8,6) = Me(6,8);  Me(8,8) = 156;   Me(8,12) = -22*L;
% Me(9,3) = Me(3,9);  Me(9,5) = Me(5,9);  Me(9,9) = 156;   Me(9,11) = 22*L;
% Me(10,4) = Me(4,10); Me(10,10) = 140*Ip/A;
% Me(11,3) = Me(3,11);    Me(11,5) = Me(5,11);    Me(11,9) = Me(9,11); Me(11,11) = 4*L^2;
% Me(12,2) = Me(2,12);    Me(12,6) = Me(6,12);    Me(12,8) = Me(8,12); Me(12,12) = 4*L^2;
% 
% Me = (rho*A*L/420)*Me;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Form the orthotropic stiffness (E1*A, E2*Iy, E3*Iz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EA  = EA;
EIy = FlapEI;  
EIz = EdgeEI;
GJ = GJ;

Ke = zeros(12,12);
Ke(1,1) = EA/L;    Ke(1,7) = -Ke(1,1);
Ke(2,2) = 12*EIz/L^3;  Ke(2,6) = 6*EIz/L^2;   Ke(2,8) = -Ke(2,2); Ke(2,12) = Ke(2,6);
Ke(3,3) = 12*EIy/L^3;  Ke(3,5) = -6*EIy/L^2;  Ke(3,9) = -Ke(3,3); Ke(3,11) = Ke(3,5);
Ke(4,4) = GJ/L;    Ke(4,10) = -Ke(4,4);
Ke(5,3) = Ke(3,5);  Ke(5,5) = 4*EIy/L; Ke(5,9) = -Ke(5,3); Ke(5,11) = 2*EIy/L;
Ke(6,2) = Ke(2,6);  Ke(6,6) = 4*EIz/L; Ke(6,8) = -Ke(6,2); Ke(6,12) = 2*EIz/L;
Ke(7,1) = -EA/L;   Ke(7,7) = EA/L;
Ke(8,2) = Ke(2,8);  Ke(8,6) = Ke(6,8);  Ke(8,8) = 12*EIz/L^3;  Ke(8,12) = Ke(8,6);
Ke(9,3) = Ke(3,9);  Ke(9,5) = Ke(5,9);  Ke(9,9) = 12*EIy/L^3;  Ke(9,11) = Ke(9,5);
Ke(10,4)= Ke(4,10); Ke(10,10) = -Ke(10,4);
Ke(11,3)= Ke(3,11); Ke(11,5) = Ke(5,11);    Ke(11,9) = Ke(9,11);    Ke(11,11) = 4*EIy/L;
Ke(12,2)= Ke(2,12); Ke(12,6) = Ke(6,12);    Ke(12,8) = Ke(8,12);    Ke(12,12) = 4*EIz/L;

% Handling Twist correction using Euler rotation sequence

rotTSmall=[1    0   0
            0   cosd(Twist) sind(Twist)
            0  -sind(Twist) cosd(Twist)];
        
Z=zeros(3);
rotT=[rotTSmall Z Z Z
        Z rotTSmall Z Z 
        Z Z rotTSmall Z
        Z Z Z rotTSmall];

Ke=rotT'*Ke*rotT;
    
end

