clc,clear

%% Constants and constraints

% init variables and symbols
S = get_variables();
syms m1 m2 m3 J3g k1 k2 k3 k4 k5 kt3 c1 c2 c3 c4 c5 ct3 L Lg theta0
var = [m1; m2; m3; J3g; k1; k2; k3; k4; k5; kt3; c1; c2; c3; c4; c5; ct3; L; Lg];

% define generalized coordinates and derivatives
syms x0 x1 x2 x3 y3 theta3 t
syms x0d x1d x2d x3d y3d theta3d 
syms x1dd x2dd x3dd y3dd theta3dd

q = [x1; x2; x3; y3; theta3];
qd = [x1d; x2d; x3d; y3d; theta3d];
qdd = [x1dd; x2dd; x3dd; y3dd; theta3dd];

% define equilibrium vectors
q_eq = [0; 0; 0; 0; 0];
qd_eq = [0; 0; 0; 0; 0;];
qdd_eq = [0; 0; 0; 0; 0;];

% define linearized coordinates and derivatives
syms x1_ x2_ x3_ y3_ theta3_
syms x1d_ x2d_ x3d_ y3d_ theta3d_
syms x1dd_ x2dd_ x3dd_ y3dd_ theta3dd_

q_ = [x1_; x2_; x3_; y3_; theta3_];
qd_ = [x1d_; x2d_; x3d_; y3d_; theta3d_];
qdd_ = [x1dd_; x2dd_; x3dd_; y3dd_; theta3dd_];

%% Constraints

% shoulder to elbow
x4 = x3 + L*cos(theta3+theta0);
x4d = x3d - L*theta3d*sin(theta3+theta0);
y4 = y3 + L*sin(theta3+theta0)-L*sin(theta0);
y4d = y3d + L*theta3d*cos(theta3+theta0);

% CoM of upper arm to elbow
x5 = x3 + Lg*cos(theta3+theta0);
x5d = x3d - Lg*theta3d*sin(theta3+theta0);
y5 = y3 + Lg*sin(theta3+theta0)-Lg*sin(theta0);
y5d = y3d + Lg*theta3d*cos(theta3+theta0);

%% Energies

% Kinetic energy
T1 = 0.5*m1*x1d^2; % wrist
T2 = 0.5*m2*x2d^2; % lower arm
T3 = 0.5*m3*x5d^2 + 0.5*m3*y5d^2 + 0.5*J3g*theta3d^2; % upper arm
T = T1 + T2 + T3;

% Potential energy
V1 = 0.5*k1*(x1 - x0)^2; % device-palm
V2 = 0.5*k2*(x2 - x1)^2; % palm-lower arm
V3 = 0.5*k3*(x3 - x2)^2 + 0.5*kt3*(theta3)^2; % lower arm-elbow
V4 = 0.5*k4*x4^2 + 0.5*k5*y4^2; % upper arm-shoulder
V = V1 + V2 + V3 + V4;

% Dissipation energy
D1 = 0.5*c1*(x1d - x0d)^2; % device-palm
D2 = 0.5*c2*(x2d - x1d)^2; % palm-lower arm
D3 = 0.5*c3*(x3d - x2d)^2 + 0.5*ct3*(theta3d)^2; % elbow
D4 = 0.5*c4*x4d^2 + 0.5*c5*y4d^2; % upper arm-shoulder
D = D1 + D2 + D3 + D4;

%% Lagrangian --> EoM

Tqd = simplify(jacobian(T,qd)).';
L1 = simplify(jacobian(Tqd, t)) + simplify( jacobian(Tqd,q)*qd + jacobian(Tqd,qd)*qdd );
L2 = simplify( jacobian(T,q) ).';
L3 = simplify( jacobian(V,q) ).';
L4 = -simplify( jacobian(D,qd) ).';

EoM = L1 - L2 + L3 - L4;

%% Linearization
M = simplify(jacobian(L1 - L2, qdd));
C = simplify(jacobian(L1 - L2 - L4, qd));
K = simplify(jacobian(EoM, q));

M_eq = subs_eq(M, q, qd, qdd, q_eq, qd_eq, qdd_eq); M_eq = subs(M_eq, theta0, pi/2);
C_eq = subs_eq(C, q, qd, qdd, q_eq, qd_eq, qdd_eq); C_eq = subs(C_eq, theta0, pi/2);
K_eq = subs_eq(K, q, qd, qdd, q_eq, qd_eq, qdd_eq); K_eq = subs(K_eq, theta0, pi/2);

Lin_EoM = M_eq*qdd_ + C_eq*qd_ + K_eq*q_;

% generalized force vector
Q = simplify(-jacobian(V, q) + -jacobian(D, qd)).';
Q = subs_eq(Q, q, qd, qdd, q_eq, qd_eq, qdd_eq);

%% Eigenmodes 

% insert variables
Kv = double(subs(K_eq,var,S));
Mv = double(subs(M_eq,var,S));
Cv = double(subs(C_eq,var,S));


% Case 1: no damping
[X,omegasquare]=eig(Kv,Mv);
omega = sqrt(diag(omegasquare));
freq = omega/(2*pi);

% Case 2: damping
[Xd,lambdad]=polyeig(Kv,Cv,Mv);
freqd = imag(lambdad)/(2*pi);

% Good approximation/ allowed if:
% 1 lightly damped
% 2 Eigenfrequencies well seperated

Beta = zeros(5,1);
Mu = zeros(5,1);
Eps = zeros(5,1);

for k = 1:size(X)
    Beta(k) = X(:,k)'*Cv*X(:,k);
    Mu(k) = transpose(X(:,k))*Mv*X(:,k);
    Eps(k) = Beta(k)/(2*omega(k)*Mu(k));
end

%% functions
function [res] = subs_eq(mat, v, vd, vdd, v_eq, vd_eq, vdd_eq)
   r1 = subs(mat, v, v_eq);
   r2 = subs(r1, vd, vd_eq);
   res = simplify( subs(r2, vdd, vdd_eq) );
end

function [S] = get_variables()
%% list of all variables
m1 = 0.45;      % kg, mass of palm
m2 = 1.15;      % kg, mass of lower arm
m3 = 1.9;       % kg, mass of upper arm
J3g = 0.0149;   % kgm^2, mass moment of inertia of upper arm
k1 = 155.8E3;   % N/m, spring constant of device-palm
k2 = 23.6E3;    % N/m, spring constant of palm-lower arm
k3 = 444.6E3;   % N/m, spring constant of lower arm- elbow
k4 = 415.4E3;   % N/m, spring constant of upper arm-shoulder, x-direction
k5 = 50.25E3;   % N/m, spring constant of upper arm-shoulder, y-direction
kt3 = 2;        % Nm/rad, rotational stiffness of elbow
c1 = 30;        % Ns/m, damping of device-palm
c2 = 202.84;    % Ns/m, damping of palm-lower arm
c3 = 500;       % Ns/m, damping of lower arm-elbow
c4 = 164.59;    % Ns/m, damping of upper arm-shoulder, x-direction
c5 = 49.99;     % Ns/m, damping of upper arm-shoulder, y-direction
ct3 = 4.9;      % Nms/rad, rotational damping of elbow
L = 0.298;      % m, length of upper arm
Lg = 0.6*L;       % m, distance between elbow and mass center of upper arm
S = [m1; m2; m3; J3g; k1; k2; k3; k4; k5; kt3; c1; c2; c3; c4; c5; ct3; L; Lg];
end