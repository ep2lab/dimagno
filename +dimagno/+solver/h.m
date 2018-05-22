function h = h(data,w,q_m,Bz,axi_r)

h = w*q_m*Bz *data.dimagno.B_on_ions + w^2*axi_r;
