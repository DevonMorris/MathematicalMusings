function [wx, wy, wz, psi, theta, phi] = euler(wx0, wy0, wz0, psi0, theta0, phi0,t,Mx,My,Mz)
  % Allocate memory
  wx = zeros(size(t));
  wy = zeros(size(t));
  wz = zeros(size(t));
  psi = zeros(size(t));
  theta = zeros(size(t));
  phi = zeros(size(t));
  dpsi = zeros(size(t));
  dtheta = zeros(size(t));
  dphi = zeros(size(t));

  % convert to radians
  psi0 = psi0*pi/180;
  theta0 = theta0*pi/180;
  phi0 = phi0*pi/180;
  wx0 = wx0*pi/180;
  wy0 = wy0*pi/180;
  wz0 = wz0*pi/180;

  % Prime the integrator
  psi(1) = psi0;
  theta(1) = theta0;
  phi(1) = phi0;
  dphi(1) = (wy0*sin(theta(1))*sin(phi(1)) + wz0*sin(theta(1))*cos(phi(1)))/cos(theta(1)) + wx0;
  dtheta(1) = wy0*cos(phi(1)) - wz0*sin(phi(1));
  dpsi(1) = (wy0*sin(phi(1)) + wz0*cos(phi(1)))/cos(theta(1));


  for i = 2:length(t)
    % Do simple RK4 integration
    dt = t(i) - t(i-1);
    q = [phi(i-1); theta(i-1); psi(i-1); dphi(i-1); dtheta(i-1); dpsi(i-1)];
    u = [Mx(i-1); My(i-1); Mz(i-1)];
    u1 = [Mx(i); My(i); Mz(i)];
    u12 = (u + u1)/2;

    k1 = dt*odefcn(t(i), q, u);
    k2 = dt*odefcn(t(i) + dt/2, q + k1/2, u12);
    k3 = dt*odefcn(t(i) + dt/2, q + k2/2, u12);
    k4 = dt*odefcn(t(i) + dt, q + k3, u1);

    q1 = q + (k1 + 2*k2 + 2*k3 + k4)/6;

    % pack up result
    phi(i) = q1(1);
    theta(i) = q1(2);
    psi(i) = q1(3);
    dphi(i) = q1(4);
    dtheta(i) = q1(5);
    dpsi(i) = q1(6);
  end

  % Convert back to angular velocities
  wx = -dpsi.*sin(theta) + dphi;
  wy = dpsi.*cos(theta).*sin(phi) + dtheta.*cos(phi);
  wz = dpsi.*cos(theta).*cos(phi) - dtheta.*sin(phi);

  % Convert back to degrees
  psi = psi*180/pi;
  theta = theta*180/pi;
  phi = phi*180/pi;
  wx = wx*180/pi;
  wy = wy*180/pi;
  wz = wz*180/pi;
end

function dqdt = odefcn(t,q,u)
  % Parameters
  I_xx = 40481.983;
  I_xy = 0;
  I_xz = 0;
  I_yy = 90353.316;
  I_yz = 0;
  I_zz = 98636.935;

  % Unpack state
  phi = q(1);
  th = q(2);
  psi = q(3);
  dphi = q(4);
  dth = q(5);
  dpsi = q(6);

  cphi = cos(phi);
  sphi = sin(phi);
  cth = cos(th);
  sth = sin(th);

  % Create mass matrix
  M = zeros(3);
  M(1,1) = I_xx;
  M(1,2) = -(I_xy*cphi - I_xz*sphi);
  M(1,3) = -(I_xx*sth + I_xz*cphi*cth + I_xy*cth*sphi);
  M(2,1) = -I_xy;
  M(2,2) = I_yy*cphi + I_yz*sphi;
  M(2,3) = I_xy*sth - I_yz*cphi*cth + I_yy*cth*sphi;
  M(3,1) = -I_xz;
  M(3,2) = -(I_yz*cphi + I_zz*sphi);
  M(3,3) = I_xz*sth + I_zz*cphi*cth - I_yz*cth*sphi;

  % Solve for dynamics
  F = zeros(3,1);
  F(1) = u(1) -((I_yz*cphi^2*cth^2 - I_yz*cth^2*sphi^2 - I_xy*cphi*cth*sth + I_xz*cth*sphi*sth - I_yy*cphi*cth^2*sphi + I_zz*cphi*cth^2*sphi)*dpsi^2 + (2*I_xz*cphi*sth - I_xx*cth + 2*I_xy*sphi*sth - I_yy*cphi^2*cth + I_zz*cphi^2*cth + I_yy*cth*sphi^2 - I_zz*cth*sphi^2 - 4*I_yz*cphi*cth*sphi)*dpsi*dth + (I_yz*sphi^2 - I_yz*cphi^2 + I_yy*cphi*sphi - I_zz*cphi*sphi)*dth^2);
  F(2) = u(2) -(I_xz*dphi^2 + (I_xx*cphi*cth - 2*I_xz*sth + I_yy*cphi*cth - I_zz*cphi*cth + 2*I_yz*cth*sphi)*dphi*dpsi + (2*I_yz*cphi - I_xx*sphi - I_yy*sphi + I_zz*sphi)*dphi*dth + (I_xz*sth^2 - I_xz*cphi^2*cth^2 - I_xx*cphi*cth*sth + I_zz*cphi*cth*sth - I_yz*cth*sphi*sth - I_xy*cphi*cth^2*sphi)*dpsi^2 + (I_xy*cth + I_xx*sphi*sth - I_yy*sphi*sth - I_zz*sphi*sth - I_xy*cphi^2*cth + I_xy*cth*sphi^2 + 2*I_xz*cphi*cth*sphi)*dpsi*dth + (I_xy*cphi*sphi - I_xz*sphi^2)*dth^2);
  F(3) = u(3) -(- I_xy*dphi^2 + (2*I_xy*sth - 2*I_yz*cphi*cth - I_xx*cth*sphi + I_yy*cth*sphi - I_zz*cth*sphi)*dphi*dpsi + (I_yy*cphi - I_xx*cphi - I_zz*cphi + 2*I_yz*sphi)*dphi*dth + (I_xy*cth^2*sphi^2 - I_xy*sth^2 + I_yz*cphi*cth*sth + I_xx*cth*sphi*sth - I_yy*cth*sphi*sth + I_xz*cphi*cth^2*sphi)*dpsi^2 + (I_xz*cth + I_xx*cphi*sth - I_yy*cphi*sth - I_zz*cphi*sth + I_xz*cphi^2*cth - I_xz*cth*sphi^2 + 2*I_xy*cphi*cth*sphi)*dpsi*dth + (I_xy*cphi^2 - I_xz*cphi*sphi)*dth^2);
  ddq = M\F;
  dqdt = zeros(6,1);

  % Pack up state derivative
  dqdt(1) = dphi;
  dqdt(2) = dth;
  dqdt(3) = dpsi;
  dqdt(4) = ddq(1);
  dqdt(5) = ddq(2);
  dqdt(6) = ddq(3);
end
