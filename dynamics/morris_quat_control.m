function [wx, wy, wz, psi, theta, phi, th1, th2, th3] = morris_quat_control(t, wx0, wy0, wz0, q0, qr)
  wx = zeros(size(t));
  wy = zeros(size(t));
  wz = zeros(size(t));
  psi = zeros(size(t));
  theta = zeros(size(t));
  phi = zeros(size(t));
  th1 = zeros(size(t));
  th2 = zeros(size(t));
  th3 = zeros(size(t));

  % Parameters
  I_xx = 40823.073;
  I_xy = -1537.807;
  I_xz = 3179.297;
  I_yy = 90593.489;
  I_yz = -128.577;
  I_zz = 98742.852;
  I = [I_xx -I_xy -I_xz; -I_xy I_yy -I_xz; -I_xz -I_yz I_zz];

  K_pq = 0.4;
  K_pw = 4;

  % store variables
  wx(1) = wx0*pi/180;
  wy(1) = wy0*pi/180;
  wz(1) = wz0*pi/180;
  q = q0;

  for i = 2:length(t)
    w = [wx(i-1); wy(i-1); wz(i-1)];
    th = [th1(i-1); th2(i-1); th3(i-1)];
    delR = logm(q.R'*qr.R);
    R = q.inv.R;
    delq = [delR(3,2); -delR(3,1); delR(2,1)]
    norm(delq);

    % Compute Control
    wr = K_pq*delq;
    delw = wr - w;
    M = cross(w, I*w) + K_pw*I*delw;

    % Just do euler integration because it's easier
    dt = t(i) - t(i-1);
    w = w + dt*dwdt(w,M)
    th = th + dt*dthdt(th,M);
    q = q*UnitQuaternion.omega(0.5*w*dt);

    % pack up states
    eul = q.torpy('deg');
    wx(i) = w(1);
    wy(i) = w(2);
    wz(i) = w(3);
    th1(i) = th(1);
    th2(i) = th(2);
    th3(i) = th(3);
    psi(i) = eul(3);
    theta(i) = eul(2);
    phi(i) = eul(1);
  end
  wx = wx*180/pi;
  wy = wy*180/pi;
  wz = wz*180/pi;
  th1 = th1*180/pi;
  th2 = th2*180/pi;
  th3 = th3*180/pi;

end

function dw = dwdt(w,M)
  % Parameters
  I_xx = 40823.073;
  I_xy = -1537.807;
  I_xz = 3179.297;
  I_yy = 90593.489;
  I_yz = -128.577;
  I_zz = 98742.852;
  I = [I_xx -I_xy -I_xz; -I_xy I_yy -I_xz; -I_xz -I_yz I_zz];
  dw = I\(M - cross(w,I*w));
end

function dth = dthdt(th, M)
  I_R = 11.9056;
  Om = 6600*2*pi/60;

  A = [-sin(th(1)) 0 cos(th(3));...
       cos(th(1)) -sin(th(2)) 0;...
       0 cos(th(2)) -sin(th(3))];

  dth = A\M;
  dth = -dth/Om;
  dth = dth/I_R;
end
