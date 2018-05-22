function Lo = charline_o(zo,ro,uo,vo)
% Returns a struct with the charline, given the data for a point of
% the charline.
% MMM20130204 
n = sqrt(uo.^2+vo.^2); % to normalize
Lo.A = uo./n;
Lo.B = zo;
Lo.C = vo./n;
Lo.D = ro;