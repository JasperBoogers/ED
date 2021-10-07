clear

% init variables and symbols
load variables.mat
syms x0 x1 x2 x3 y3 theta3 t
syms x0d x1d x2d x3d y3d theta3d 
syms x1dd x2dd x3dd y3dd theta3dd

% define generalized coordinates vectors
q = [x1; x2; x3; y3; theta3];
qd = [x1d; x2d; x3d; y3d; theta3d];
qdd = [x1dd; x2dd; x3dd; y3dd; theta3dd];

% define constraints and corresponding derivatives
x4 = x3 + L*cos(theta3);
x4d = x3d - L*theta3d*sin(theta3);
y4 = y3 + L*sin(theta3);
y4d = y3d + L*theta3d*cos(theta3);

% define kinetic energy
T1 = 0.5*m1*x1d^2; % wrist
T2 = 0.5*m2*x2d^2; % lower arm
T3 = 0.5*m3*x3d^2 + 0.5*m3*y3d^2 + 0.5*(J3g + m3*Lg^2)*theta3d^2; % upper arm
T = T1 + T2 + T3;

% define potential energy
V1 = 0.5*k1*(x1 - x0)^2; % device-palm
V2 = 0.5*k2*(x2 - x1)^2; % palm-lower arm
V3 = 0.5*k3*((x3 - x2)^2 + y3^2) + 0.5*kt3*theta3^2; % lower arm-elbow
V4 = 0.5*k4*(x4^2 + y4^2) + 0.5*k5*(x4^2 + y4^2); % upper arm-shoulder
V = V1 + V2 + V3 + V4;

% define dissipation
D1 = 0.5*c1*(x1d - x0d)^2; % device-palm
D2 = 0.5*c2*(x2d - x1d)^2; % palm-lower arm
D3 = 0.5*c3*( (x3d - x2d)^2 + y3d^2) + 0.5*ct3*theta3d^2; % elbow
D4 = 0.5*c4*x4d^2 + 0.5*c5*y4d^2; % upper arm-shoulder
D = D1 + D2 + D3 + D4;

% construct lagrangian
Tqd = simplify(jacobian(T,qd))';
L1 = simplify(jacobian(Tqd, t)) + simplify( jacobian(Tqd,q)*qd + jacobian(Tqd,qd) )*qdd;
L2 = simplify( jacobian(T,q) ).';
L3 = simplify( jacobian(V,q) ).';
L4 = simplify( jacobian(D,qd) ).';

EoM = L1 - L2 + L3 - L4;