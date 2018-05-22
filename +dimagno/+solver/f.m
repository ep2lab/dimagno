function f = f(data,u,v,Bz,Br,ne,we,axi_r)

f = 0;
temp = (Bz*v - Br*u);
for i = 1:data.plasma.n_electrons
    qe = data.plasma.electrons{i}.q;    me = data.plasma.electrons{i}.m;
    f = f - (qe*we{i}*temp - me*we{i}^2*axi_r*v)...
        /data.plasma.electrons{i}.dh_dn(ne{i});
end
f = f/sum(cell2mat(ne)) - v*axi_r;
    
    
    
