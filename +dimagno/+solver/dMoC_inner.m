%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p4,iter,iter_err] = dMoC_inner(data,i1,i2,F)

% abbrv
m = data.plasma.ions.m;
q = data.plasma.ions.q;
q_m = q/m;

%% Get p1, p2

p1 = F.point(i1);
p2 = F.point(i2);
% Complete p1, p2
p1.Hi = q*p1.phi + 0.5 *m* (p1.u.^2 + p1.v.^2 + p1.w.^2);
p2.Hi = q*p2.phi + 0.5 *m* (p2.u.^2 + p2.v.^2 + p2.w.^2); 
 
% Variables useful for determining 3 (see NBII p105)
p3.lzd_ = p2.u-p1.u; p3.lrd_ = p2.v-p1.v;
p3.zd_ = p2.z-p1.z; p3.rd_ = p2.r-p1.r;
p3.a_ = p3.lzd_*p3.rd_-p3.lrd_*p3.zd_;
 
%% Main iteration routine

if data.dimagno.corrector % (for corrector, if is_corrector_on == 1)
    max_iter = data.dimagno.max_iter;
    Z = Inf; R = Inf; U = Inf;V = Inf;W = Inf;PHI = Inf; % convergence checkers
else
    max_iter = 0;
    iter_err = NaN;
end

for iter = 0:max_iter
    
    %% Properties for p,m
    
    if iter == 0 
        pp = p2;
        pm = p1;
    else
        pp.z = (p2.z+p4.z)/2; pp.r = (p2.r+p4.r)/2;
        pm.z = (p1.z+p4.z)/2; pm.r = (p1.r+p4.r)/2;
                 
        pp.u = (p2.u+p4.u)/2; pp.v = (p2.v+p4.v)/2;
        pm.u = (p1.u+p4.u)/2; pm.v = (p1.v+p4.v)/2;
    
        pp.phi = (p2.phi+p4.phi)/2;
        pm.phi = (p1.phi+p4.phi)/2;
      
        [pp.psi,pp.Bz,pp.Br] = data.field.field_2d(pp.z,pp.r);
        [pm.psi,pm.Bz,pm.Br] = data.field.field_2d(pm.z,pm.r);
        pp.w = (p2.w+p4.w)/2;
        pm.w = (p1.w+p4.w)/2;
        pp.etae = data.ic.calculate_etae_2d(pp.r,pp.psi);
        pm.etae = data.ic.calculate_etae_2d(pm.r,pm.psi);
        pp.We = data.ic.get_We(pp.etae);
        pm.We = data.ic.get_We(pm.etae);
        pp.He = data.ic.get_He(pp.etae);
        pm.He = data.ic.get_He(pm.etae);
        for j = 1:data.plasma.n_electrons
            pp.we{j} = pp.We{j}*pp.r;
            pm.we{j} = pm.We{j}*pm.r;
            pp.ne{j} = data.plasma.electrons{j}.n(pp.He{j}-data.plasma.electrons{j}.q*pp.phi);
            pm.ne{j} = data.plasma.electrons{j}.n(pm.He{j}-data.plasma.electrons{j}.q*pm.phi);
        end        
        pp.c = data.plasma.cs(pp.ne);
        pm.c = data.plasma.cs(pm.ne);
    end
    
    % Slope terms for lines p,m
    
    Lp = dimagno.solver.charline_p(p2.z,p2.r,pp.u,pp.v,pp.c);
    Lm = dimagno.solver.charline_m(p1.z,p1.r,pm.u,pm.v,pm.c);
     
    %% Obtain position of point 4 and line parameters tm, tp

    denom = Lp.A.*Lm.C - Lp.C.*Lm.A;
    tm = (Lp.C.*(Lm.B-Lp.B) - Lp.A.*(Lm.D-Lp.D)) ./ denom;
    tp = - (Lm.C.*(Lp.B-Lm.B) - Lm.A.*(Lp.D-Lm.D)) ./ denom;
    p4.z = Lm.A*tm+Lm.B;
    p4.r = Lm.C*tm+Lm.D;
    [p4.psi,p4.Bz,p4.Br] = data.field.field_2d(p4.z,p4.r);
    
    %% Obtain position of point 3
    
    % Variables needed to find p3 (see NBII p105)
    z41_ = p4.z-p1.z;
    r41_ = p4.r-p1.r;
    p3.b_ = p3.lrd_*z41_ - p3.lzd_*r41_ - p1.v*p3.zd_ + p1.u*p3.rd_;
    p3.c_ = p1.v*z41_ - p1.u*r41_;
    s3 = -2*p3.c_/(p3.b_*(1+sqrt(1-4*p3.a_*p3.c_/p3.b_^2)));
    
    % Interpolate for properties at 3 needed for slope terms
    p3.z = p1.z*(1-s3) + p2.z*s3;
    p3.r = p1.r*(1-s3) + p2.r*s3;
    p3.u = p1.u*(1-s3) + p2.u*s3;
    p3.v = p1.v*(1-s3) + p2.v*s3;
    p3.w = p1.w*(1-s3) + p2.w*s3;
    [p3.psi,~,~] = data.field.field_2d(p3.z,p3.r);
    p3.Hi = p1.Hi*(1-s3) + p2.Hi*s3;
    p3.phi = p3.Hi/q - 0.5 * (p3.u.^2 + p3.v.^2 + p3.w.^2)/q_m;
      
    % Properties at po
    if iter == 0
        po = p3;
    else
        po.u = (p3.u + p4.u)/2; po.v = (p3.v + p4.v)/2;        
    end
    
    %% canonical momentum equation on charline o (needs nothing more) --> p4.w
    
    if data.field.is_axi
        p4.w = (p3.r*p3.w + q_m*(p3.psi-p4.psi)*data.dimagno.B_on_ions )/p4.r;
    else
        p4.w = p3.w + q_m*(p3.psi-p4.psi)*data.dimagno.B_on_ions;
    end
    
    %% Solve equations to calculate properties at 4
    
    % f,g,h forcing terms on each characteristic line
    axi_rp = 0;
    if data.field.is_axi && pp.r ~= 0
        axi_rp = 1/pp.r; % This adds the centrifugal contribution iff axi == 1 and r_ ~= 0
    end
    fp = dimagno.solver.f(data,pp.u,pp.v,pp.Bz,pp.Br,pp.ne,pp.we,axi_rp);
    gp = dimagno.solver.g(data,pp.w,q_m,pp.Br);
    hp = dimagno.solver.h(data,pp.w,q_m,pp.Bz,axi_rp);
    
    axi_rm = 0;
    if data.field.is_axi
        axi_rm = 1/pm.r;
    end
    fm = dimagno.solver.f(data,pm.u,pm.v,pm.Bz,pm.Br,pm.ne,pm.we,axi_rm);
    gm = dimagno.solver.g(data,pm.w,q_m,pm.Br);
    hm = dimagno.solver.h(data,pm.w,q_m,pm.Bz,axi_rm);
  
    % System Matrix
    MC = [pp.v, -pp.u, -q_m*sqrt((pp.u^2+pp.v^2)/pp.c^2-1); % char. p
          pm.v, -pm.u, +q_m*sqrt((pm.u^2+pm.v^2)/pm.c^2-1); % char. m
          po.u,  po.v,                                 q_m]; % char. o
    
    % Column vector containing values at previous points
    Gi(1,1) = MC(1,:)*[p2.u;p2.v;p2.phi];
    Gi(2,1) = MC(2,:)*[p1.u;p1.v;p1.phi];
    Gi(3,1) = MC(3,:)*[p3.u;p3.v;p3.phi];
    
    % Column vector containing forcing terms (multiplied by the line parameter increments from points 2,1,3, respectively
    Fi(1,1) = ((Lp.A*pp.v-Lp.C*pp.u)*fp + Lp.C*gp - Lp.A*hp) * tp;
    Fi(2,1) = ((Lm.A*pm.v-Lm.C*pm.u)*fm + Lm.C*gm - Lm.A*hm) * tm;
    Fi(3,1) = 0.5*(p3.w^2-p4.w^2);
     
    % Solve for properties at point 4
    sol = MC\(Gi+Fi);
    p4.u   = sol(1); 
    p4.v   = sol(2);
    
    %% Enforce exact energy equation along char o (after u,v at 4 are known)
    
    p4.phi = p3.Hi/q - 0.5*(p4.u^2+p4.v^2+p4.w^2)/q_m;
    
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
    error([mfilename,':NOCONV'],'Convergence was not achieved for p4:[z,r] = [%G,%G]',p4.z,p4.r);
end

%% Additional properties at 4

p4.etae = data.ic.calculate_etae_2d(p4.r,p4.psi);
p4.He = data.ic.get_He(p4.etae);
p4.We = data.ic.get_We(p4.etae);
for i = 1:data.plasma.n_electrons
    p4.we{i} = p4.We{i}*p4.r;
    p4.ne{i} = data.plasma.electrons{i}.n(p4.He{i}-data.plasma.electrons{i}.q*p4.phi);
    p4.Te{i} = data.plasma.electrons{i}.T(p4.ne{i});
end
p4.c = data.plasma.cs(p4.ne);
