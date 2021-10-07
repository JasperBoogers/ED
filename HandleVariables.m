clear

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


%% save as .mat file
save -mat variables.mat