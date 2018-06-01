%{
Function to calculate variables at a new inner point using
direct MOC (p1 and p2 are given; p3 is found and interpolated).

input:
******

data: data structure
i1,i2,F: indices of p1 and p2 in front F

output:
*******

p4 (structure): the new point with all properties (basic and derived).
iter,iter_err: convergence information

MMM20130310
%}
function [p4,iter,iter_err] = dMoC_axis(data,i1,i3,F)

% abbrv
m = data.plasma.ions.m;
q = data.plasma.ions.q;
q_m = q/m;

%% Get p1, p3

p1 = F.point(i1);
p3 = F.point(i3);
% Ion energy for point 3 calculations, so that it is better conserved 20120424
p3.Hi = q*p3.phi + 0.5 *m* (p3.u.^2 + p3.v.^2 + p3.w.^2);

% Char line o
Lo = dimagno.straightline;
Lo.A = p3.z.*0+1; Lo.B = p3.z; Lo.C = p3.z*0; Lo.D = p3.r;

%% Already known p4 properties

p4.r = 0;
p4.v = 0;
p4.w = 0;
for i = 1:data.plasma.n_electrons
    p4.we{i} = 0;
end
p4.psi = 0;
p4.Br = 0;
 
%% Main iteration routine

if data.dimagno.corrector % (for corrector, if is_corrector_on == 1)
    max_iter = data.dimagno.max_iter;
    Z = Inf; R = Inf; U = Inf;V = Inf;W = Inf;PHI = Inf; % convergence checkers
else
    max_iter = 0;
    iter_err = NaN;
end

for iter = 0:max_iter
    
    %% Properties for o,m
    
    if iter == 0 
        pm = p1;
    else
        pm.z = (p1.z+p4.z)/2; pm.r = (p1.r+p4.r)/2;
        pm.u = (p1.u+p4.u)/2; pm.v = (p1.v+p4.v)/2;
        pm.phi = (p1.phi+p4.phi)/2;
      
        [pm.psi,pm.Bz,pm.Br] = data.field.field_2d(pm.z,pm.r);
        pm.w = (p1.w+p4.w)/2;
        pm.etae = data.ic.calculate_etae_2d(pm.r,pm.psi);
        pm.We = data.ic.get_We(pm.etae);
        pm.He = data.ic.get_He(pm.etae);
        for j = 1:data.plasma.n_electrons
            pm.we{j} = pm.We{j}*pm.r;
            pm.ne{j} = data.plasma.electrons{j}.n(pm.He{j}-data.plasma.electrons{j}.q*pm.phi);
        end        
        pm.c = data.plasma.cs(pm.ne);
    end
    
    % Slope terms for lines p,m
    
    Lmm = dimagno.solver.charline_m(p1.z,p1.r,pm.u,pm.v,pm.c);
    Lm = dimagno.straightline(Lmm);
     
    % Properties at po
    if iter == 0
        po = p3;
    else
        po.u = (p3.u + p4.u)/2; po.v = (p3.v + p4.v)/2;        
    end
    
    %% Obtain position of point 4 and line parameters tm, tp
    
    [p4.z,p4.r,tm] = Lm.intersectionxy(Lo);
    [~,p4.Bz] = data.field.field_2d(p4.z,p4.r);
    
    %% Solve equations to calculate properties at 4
    
    % f,g,h forcing terms on each characteristic line
    axi_rm = 0;
    if data.field.is_axi
        axi_rm = 1/pm.r;
    end
    fm = dimagno.solver.f(data,pm.u,pm.v,pm.Bz,pm.Br,pm.ne,pm.we,axi_rm);
    gm = dimagno.solver.g(data,pm.w,q_m,pm.Br);
    hm = dimagno.solver.h(data,pm.w,q_m,pm.Bz,axi_rm);
         
    % System Matrix
    MC = [pm.v, +q_m*sqrt(pm.u^2+pm.v^2-pm.c^2)/pm.c; % char. m
          po.u,                                 q_m]; % char. o
    
    % Column vector containing values at previous points
    Gi(1,1) = MC(1,:)*[p1.u;p1.phi] - pm.u*p1.v; % this is the contribution of the missing column in the matrix MC
    Gi(2,1) = MC(2,:)*[p3.u;p3.phi];
    
    % Column vector containing forcing terms (multiplied by the line parameter increments from points 2,1,3, respectively
    Fi(1,1) = ((Lm.A*pm.v-Lm.C*pm.u)*fm + Lm.C*gm - Lm.A*hm) * tm;
    Fi(2,1) = 0;
    
    % Solve for properties at point 4
    sol = MC\(Gi+Fi);
    p4.u   = sol(1);
    
    %% Enforce exact energy equation along char o (after u,v at 4 are known)
    
    p4.phi = p3.Hi/q - 0.5*p4.u^2/q_m;
    
    %% Convergence check (if corrector is on)
    
    if ~data.dimagno.corrector
        break; % no corrector steps
    else
        iter_err = abs(p4.z-Z)+abs(p4.r-R)+abs(p4.u-U)+abs(p4.v-V)+abs(p4.w-W)+abs(p4.phi-PHI);
        if ( iter_err < data.dimagno.tol )
            break;
        else
            Z = p4.z;R = p4.r;U = p4.u;V = p4.v;W = p4.w;PHI = p4.phi;
        end
    end
    
end

%% If program reaches here, corrector is on and convergence was not achieved

if iter == max_iter && data.dimagno.corrector 
    error([mfilename,':NOCONV'],'Axis point convergence was not achieved for p4:[z,r] = [%G,%G]',p4.z,p4.r);
end

%% Additional properties at 4

p4.etae = data.ic.calculate_etae_2d(p4.r,p4.psi);
p4.He = data.ic.get_He(p4.etae); 
for i = 1:data.plasma.n_electrons 
    p4.ne{i} = data.plasma.electrons{i}.n(p4.He{i}-data.plasma.electrons{i}.q*p4.phi);
    p4.Te{i} = data.plasma.electrons{i}.T(p4.ne{i});
end
p4.c = data.plasma.cs(p4.ne);
