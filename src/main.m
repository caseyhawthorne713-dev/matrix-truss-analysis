clc; clear; close all;

% Define nodes, elements, material, and loads
L = 1;                  % Length unit in feet
A = 4/144;              % Bar cross-sectional area (ft^2)
E = 30e6;               % Modulus of elasticity (psi)
Q = 1000;               % Applied load (lbs)

% Node coordinates (x, y)
Coords = [ 
    0 0;      %1
    3*L 0;    %2
    27*L 0;   %3
    30*L 0;   %4
    3*L 4*L;  %5
    6*L 4*L;  %6
    9*L 4*L;  %7
    12*L 4*L; %8
    15*L 4*L; %9
    18*L 4*L; %10
    21*L 4*L; %11
    24*L 4*L; %12
    27*L 4*L; %13
    6*L 8*L;  %14
    12*L 8*L; %15
    18*L 8*L; %16
    24*L 8*L; %17
];

% Applied nodal forces (Fx, Fy)
F_Applied = [ 
    0 0;    %1
    0 0;    %2
    0 0;    %3
    0 0;    %4
    0 0;    %5
    0 0;    %6
    0 0;    %7
    0 -Q;   %8
    0 0;    %9
    0 -2*Q; %10
    0 0;    %11
    0 0;    %12
    0 0;    %13
    0 0;    %14
    0 0;    %15
    0 0;    %16
    0 0;    %17
];

% Boundary conditions (1=fixed, 0=free)
Bounds = [ 
    1 1; %1
    1 1; %2
    1 1; %3
    1 1; %4
    0 0; %5
    0 0; %6
    0 0; %7
    0 0; %8
    0 0; %9
    0 0; %10
    0 0; %11
    0 0; %12
    0 0; %13
    0 0; %14
    0 0; %15
    0 0; %16
    0 0; %17
];

% Prescribed displacements (0 for all)
U_prescribed = zeros(size(Bounds));

% Element connectivity (node pairs)
Connectivity = [
    1 5;
    2 5;
    2 6;
    12 3;
    13 3;
    13 4;
    5 6;
    6 7;
    7 8;
    8 9;
    9 10;
    10 11;
    11 12;
    12 13;
    5 14;
    14 6;
    14 7;
    7 15;
    15 8;
    15 9;
    9 16;
    16 10;
    16 11;
    11 17;
    17 12;
    17 13;
    14 15;
    15 16;
    16 17;
];

nElem = size(Connectivity,1);
elem_A = A*ones(nElem,1);
elem_E = E*ones(nElem,1);

nNodes = size(Coords,1);
GDof = 2*nNodes; % 2 DOFs per node

% Global matrices
K = zeros(GDof, GDof);                 % Global stiffness
F = reshape(F_Applied', [], 1);        % Global force vector
U_presc = reshape(U_prescribed', [], 1);

% Identify fixed and free DOFs
fixedDofs = [];
for i = 1:nNodes
    if Bounds(i,1) == 1, fixedDofs(end+1) = 2*i-1; end
    if Bounds(i,2) == 1, fixedDofs(end+1) = 2*i; end
end
freeDofs = setdiff(1:GDof, fixedDofs);

%% Assemble global stiffness matrix
for e = 1:nElem
    n1 = Connectivity(e,1); n2 = Connectivity(e,2);
    x1 = Coords(n1,1); y1 = Coords(n1,2);
    x2 = Coords(n2,1); y2 = Coords(n2,2);
    
    dx = x2 - x1; dy = y2 - y1;
    Le = sqrt(dx^2 + dy^2);
    c = dx/Le; s = dy/Le;
    
    Ae = elem_A(e); Ee = elem_E(e);
    
    ke = (Ae*Ee/Le) * [ c^2  c*s -c^2 -c*s;
                        c*s  s^2 -c*s -s^2;
                       -c^2 -c*s  c^2  c*s;
                       -c*s -s^2  c*s  s^2 ];
    
    dof = [2*n1-1 2*n1 2*n2-1 2*n2];
    K(dof,dof) = K(dof,dof) + ke;
end

% Solve for U
U = zeros(GDof,1);
U(fixedDofs) = U_presc(fixedDofs);

Kff = K(freeDofs, freeDofs);
Kfc = K(freeDofs, fixedDofs);
Ff  = F(freeDofs);

U(freeDofs) = Kff \ (Ff - Kfc * U(fixedDofs));

% Solve for RxN at 1 2 3 4
Rc = K(fixedDofs,:) * U - F(fixedDofs);

% Compute member stresses
memberStress = zeros(nElem,1);
for e = 1:nElem
    n1 = Connectivity(e,1); n2 = Connectivity(e,2);
    dx = Coords(n2,1) - Coords(n1,1);
    dy = Coords(n2,2) - Coords(n1,2);
    Le = sqrt(dx^2 + dy^2);
    c = dx/Le; s = dy/Le;
    Ae = elem_A(e); Ee = elem_E(e);
    dof = [2*n1-1 2*n1 2*n2-1 2*n2];
    Ue = U(dof);
    memberStress(e) = (Ee/Le) * [-c -s c s] * Ue;  % Axial stress = Force / Area
end

% Plot truss
X = Coords(:,1); Y = Coords(:,2);
scale = 5; % magnification factor
Xd = X + scale*U(1:2:end);
Yd = Y + scale*U(2:2:end);

figure; hold on; grid on;
title('Truss Deformation');
xlabel('X (ft)'); ylabel('Y (ft)');

% Plot undeformed truss
for e = 1:nElem
    n1 = Connectivity(e,1); n2 = Connectivity(e,2);
    plot([X(n1), X(n2)], [Y(n1), Y(n2)], 'b--');
end

% Plot deformed truss
for e = 1:nElem
    n1 = Connectivity(e,1); n2 = Connectivity(e,2);
    plot([Xd(n1), Xd(n2)], [Yd(n1), Yd(n2)], 'r-', 'LineWidth',2);
end

legend('Undeformed','Deformed'); axis equal;

% Results
disp('Nodal Displacements');
disp(reshape(U,2,[])');

disp('Reaction Forces at Supports');
disp(Rc);

disp('Member Axial Stresses (psi, + tension / - compression)');
disp(memberStress);
