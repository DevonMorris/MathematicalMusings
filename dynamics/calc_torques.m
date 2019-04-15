function M = calc_torques(q)
  % Parameters
  I_xx = 40823.073;
  I_xy = -1537.807;
  I_xz = 3179.297;
  I_yy = 90593.489;
  I_yz = -128.577;
  I_zz = 98742.852;

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

  % Solve for dynamics
  M = zeros(3,1);
  M(1) = ((I_yz*cphi^2*cth^2 - I_yz*cth^2*sphi^2 - I_xy*cphi*cth*sth + I_xz*cth*sphi*sth - I_yy*cphi*cth^2*sphi + I_zz*cphi*cth^2*sphi)*dpsi^2 + (2*I_xz*cphi*sth - I_xx*cth + 2*I_xy*sphi*sth - I_yy*cphi^2*cth + I_zz*cphi^2*cth + I_yy*cth*sphi^2 - I_zz*cth*sphi^2 - 4*I_yz*cphi*cth*sphi)*dpsi*dth + (I_yz*sphi^2 - I_yz*cphi^2 + I_yy*cphi*sphi - I_zz*cphi*sphi)*dth^2);
  M(2) = (I_xz*dphi^2 + (I_xx*cphi*cth - 2*I_xz*sth + I_yy*cphi*cth - I_zz*cphi*cth + 2*I_yz*cth*sphi)*dphi*dpsi + (2*I_yz*cphi - I_xx*sphi - I_yy*sphi + I_zz*sphi)*dphi*dth + (I_xz*sth^2 - I_xz*cphi^2*cth^2 - I_xx*cphi*cth*sth + I_zz*cphi*cth*sth - I_yz*cth*sphi*sth - I_xy*cphi*cth^2*sphi)*dpsi^2 + (I_xy*cth + I_xx*sphi*sth - I_yy*sphi*sth - I_zz*sphi*sth - I_xy*cphi^2*cth + I_xy*cth*sphi^2 + 2*I_xz*cphi*cth*sphi)*dpsi*dth + (I_xy*cphi*sphi - I_xz*sphi^2)*dth^2);
  M(3) = (- I_xy*dphi^2 + (2*I_xy*sth - 2*I_yz*cphi*cth - I_xx*cth*sphi + I_yy*cth*sphi - I_zz*cth*sphi)*dphi*dpsi + (I_yy*cphi - I_xx*cphi - I_zz*cphi + 2*I_yz*sphi)*dphi*dth + (I_xy*cth^2*sphi^2 - I_xy*sth^2 + I_yz*cphi*cth*sth + I_xx*cth*sphi*sth - I_yy*cth*sphi*sth + I_xz*cphi*cth^2*sphi)*dpsi^2 + (I_xz*cth + I_xx*cphi*sth - I_yy*cphi*sth - I_zz*cphi*sth + I_xz*cphi^2*cth - I_xz*cth*sphi^2 + 2*I_xy*cphi*cth*sphi)*dpsi*dth + (I_xy*cphi^2 - I_xz*cphi*sphi)*dth^2);
end
