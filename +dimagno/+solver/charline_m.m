function Lm = charline_m(zm,rm,um,vm,cm)
% Returns a struct object with the charline, given the data for a point of
% the charline.
% 20130206: I have updated the expressions according to NBII p.81 !!
% MMM20130204 
sqrtm = sqrt(um.^2+vm.^2-cm.^2);
Lm.A = um./cm+vm./sqrtm;
Lm.B = zm;
Lm.C = vm./cm-um./sqrtm;
Lm.D = rm;
% normalize
n = sqrt(Lm.A.^2 + Lm.C.^2);
Lm.A = Lm.A./n;
Lm.C = Lm.C./n;