clear all;
clc;

% Time (s)
dt = 0.01;      % Step size
tf = 30;        % Final time
t = 0:dt:tf;    % Time

% Initial conditions (deg and deg/s)
wx0 = 0;
wy0 = 0;
wz0 = 0;
q0 = Quaternion(randn(), randn(1,3));
q0 = q0*sign(q0.s);
q0 = q0.unit();
% This makes sure that q0 is on the first cover of SO(3)

% This is our reference quaternion, we will start with a random static reference
qr = Quaternion(randn(), randn(1,3));
%qr = Quaternion(1, [0 0 0]);
qr = qr*sign(qr.s);
qr = qr.unit();

% Call the project function. Note that the initial values are passed in
% units of degrees and degrees/s, and the function returns solution
% vectors in units of degrees and degrees/s.
%
% ** CHANGE THIS TO THE NAME OF YOUR OWN FUNCTION. THE FUNCTION
% PARAMETERS AND ORDER SHOULD STAY THE SAME. THE NAME OF YOUR FUNCTION
% SHOULD BE lastname1_lastname2_cmg() ** Your function must return
% values of wx, wy, wz, psi, theta, and phi at each time step defined
% in the vector t.
[wx,wy,wz,psi,theta,phi,th1, th2,th3]=morris_quat_control(t,wx0,wy0,wz0,q0,qr);

% convert qr to euler for plotting
eul = qr.torpy('deg');

% Plot results
subplot(3,1,1);
plot(t,wx,t,wy,t,wz);
xlabel('t (s)');
ylabel('\omega (deg/s)');
legend('\omega_x','\omega_y','\omega_z');
subplot(3,1,2);
plot(t,psi,t,theta,t,phi)
hold on
plot(t,ones(size(t))*eul(3), '--',t,ones(size(t))*eul(2), '--',t,ones(size(t))*eul(1), '--');
xlabel('t (s)');
ylabel('\psi, \theta, \phi (deg)');
legend('\psi','\theta','\phi', '\psi_r', '\theta_r', '\phi_r');
subplot(3,1,3);
plot(t,th1,t,th2,t,th3);
xlabel('t (s)');
ylabel('\theta_1, \theta_2, \theta_3 (deg)');
legend('\theta_1', '\theta_2', '\theta_3');
