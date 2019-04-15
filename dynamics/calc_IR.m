function IR = calc_IR(dw)
  % Parameters
  I_xx = 40823.073;
  I_xy = -1537.807;
  I_xz = 3179.297;
  I_yy = 90593.489;
  I_yz = -128.577;
  I_zz = 98742.852;

  dth = 30*pi/180;
  I = [I_xx -I_xy -I_xz; -I_xy I_yy -I_yz; -I_xy -I_yz I_zz];
  Om = 6600*2*pi/60;
  IR_vec = (I*dw)/(Om*dth);
  IR = IR_vec;
end

