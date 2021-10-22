clear

%% Constants and constraints

% init variables and symbols
S = struct2cell(load('variables.mat'));
syms m1 m2 m3 J3g k1 k2 k3 k4 k5 kt3 c1 c2 c3 c4 c5 ct3 L Lg
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

% define constraints and corresponding derivatives
x4 = x3 + L*cos(theta3);
x4d = x3d - L*theta3d*sin(theta3);
y4 = y3 + L*sin(theta3);
y4d = y3d + L*theta3d*cos(theta3);

%% Energies

% Kinetic energy
T1 = 0.5*m1*x1d^2; % wrist
T2 = 0.5*m2*x2d^2; % lower arm
T3 = 0.5*m3*x3d^2 + 0.5*m3*y3d^2 + 0.5*(J3g + m3*Lg^2)*theta3d^2; % upper arm
T = T1 + T2 + T3;

% Potential energy
V1 = 0.5*k1*(x1 - x0)^2; % device-palm
V2 = 0.5*k2*(x2 - x1)^2; % palm-lower arm
V3 = 0.5*k3*(x3 - x2)^2 + 0.5*kt3*(theta3-pi/2)^2; % lower arm-elbow
V4 = 0.5*k4*(x4^2 + (y4-L)^2) + 0.5*k5*(x4^2 + (y4-L)^2); % upper arm-shoulder
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

M_eq = subs_eq(M, q, qd, qdd, q_eq, qd_eq, qdd_eq);
C_eq = subs_eq(C, q, qd, qdd, q_eq, qd_eq, qdd_eq);
K_eq = subs_eq(K, q, qd, qdd, q_eq, qd_eq, qdd_eq);

Lin_EoM = M_eq*qdd_ + C_eq*qd_ + K_eq*q_;

% possible solution for generalized force vector
Q = simplify(-jacobian(V, q) + -jacobian(D, qd)).';
Q = subs_eq(Q, q, qd, qdd, q_eq, qd_eq, qdd_eq);

%% Eigenmodes 
% No damping
Kv = double(subs(K_eq,var,S));
Mv = double(subs(M_eq,var,S));
Cv = double(subs(C_eq,var,S));
[X,eigenfreq]=eig(Kv,Mv);
omega = sqrt(diag(eigenfreq));

% Damping
% Good approximation/ allowed if:
% 1 lightly damped
% 2 Eigenfrequencies well seperated
Beta = X'*Cv.*X;
for k = 1:size(X)
    Gamma(k) = transpose(X(:,k))*Kv*X(:,k);
    Mu(k) = transpose(X(:,k))*Mv*X(:,k);
    lambdaD = -0.5*diag(Beta)/Mu(k) + 1i*omega;
end

alpha = zeros(5);
for k = 1:size(X)
    for s=1:size(X)
        if (k~=s)
            alpha(:,k) = alpha(:,k) + 1i.*omega(k).*Beta(k,s).*X(:,s)/(Mu(s).*(omega(k)^2-omega(s)^2));
        else
           % do nothing
        end
    end
end
Z = X + alpha;


%% functions
function [res] = subs_eq(mat, v, vd, vdd, v_eq, vd_eq, vdd_eq)
   r1 = subs(mat, v, v_eq);
   r2 = subs(r1, vd, vd_eq);
   res = simplify( subs(r2, vdd, vdd_eq) );
end