%   Non-uniform, Orthotropic Beam for FAST Blade File

function [M,K] = ortho_beam_fem(nelems,props)

% Based on code by D. Todd Griffith heavily modified by Mayank Chetan
% %
% M and K are global system matrices

L   = props(:,1);
MassDen  = props(:,2);
FlapEI  = props(:,3);
EdgeEI  = props(:,4);
GJ   = props(:,5);
EA  = props(:,6);
Twist  = props(:,7);

%  Begin forming system matrices

nnodes = nelems+1;
M = zeros(6*nnodes,6*nnodes);
K = zeros(6*nnodes,6*nnodes);

[Me,Ke]=ortho_beam_elem(L(1),MassDen(1),FlapEI(1),EdgeEI(1),GJ(1),EA(1),Twist(1));
M(1:12,1:12) = Me;
K(1:12,1:12) = Ke;

ind1 = 1;

if nelems > 1
    for i = 2:nelems

        ind1 = ind1 + 6;
        ind2 = ind1 + 11;
        
        [Me,Ke]=ortho_beam_elem(L(i),MassDen(i),FlapEI(i),EdgeEI(i),GJ(i),EA(i),Twist(i));

        M(ind1:ind2,ind1:ind2) = M(ind1:ind2,ind1:ind2) + Me;
        K(ind1:ind2,ind1:ind2) = K(ind1:ind2,ind1:ind2) + Ke;
        
    end
end


end