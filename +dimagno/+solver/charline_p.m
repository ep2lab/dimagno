function Lp = charline_p(zp,rp,up,vp,cp)
% Returns a struct with the charline, given the data for a point of
% the charline.
% 20130206: I have updated the expressions according to NBII p.81 !!
% MMM20130204 
sqrtp = sqrt(up.^2+vp.^2-cp.^2);
Lp.A = up./cp-vp./sqrtp;
Lp.B = zp;
Lp.C = vp./cp+up./sqrtp;
Lp.D = rp;
% normalize
n = sqrt(Lp.A.^2 + Lp.C.^2);
Lp.A = Lp.A./n;
Lp.C = Lp.C./n;