%----------------------------------------------------------------------
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
%----------------------------------------------------------------------

function [p4,iter,iter_err] = dMoC_border(data,i3,i2,F)

% abbrv
m = data.plasma.ions.m;
q = data.plasma.ions.q;
q_m = q/m;

%% Get p1, p2

p3 = F.point(i3);
p2 = F.point(i2);
% Ion energy for point 3 calculations, so that it is better conserved 20120424
p3.Hi = q*p3.phi + 0.5 *m* (p3.u.^2 + p3.v.^2 + p3.w.^2); 
 
% Necesary variables for determining p4 position
for j = 1:data.plasma.n_electrons
    We_(j) = data.ic.We_{j}(end);
    De_(j) = data.ic.De_{j}(end);
    me_(j) = data.plasma.electrons{j}.m;
end    
[~,j_] = max(me_.*We_); % The electron species that has max slope

%% Main iteration routine

if data.dimagno.corrector % (for corrector, if is_corrector_on == 1)
    max_iter = data.dimagno.max_iter;
    Z = Inf; R = Inf; U = Inf;V = Inf;W = Inf;PHI = Inf; % convergence checkers
else
    max_iter = 0;
    iter_err = NaN;
end

for iter = 0:max_iter
    
    %% Properties for p,o    
    if iter == 0 
        pp = p2; 
    else
        pp.z = (p2.z+p4.z)/2; pp.r = (p2.r+p4.r)/2;  
        pp.u = (p2.u+p4.u)/2; pp.v = (p2.v+p4.v)/2;  
        pp.phi = (p2.phi+p4.phi)/2; 
      
        [pp.psi,pp.Bz,pp.Br] = data.field.field_2d(pp.z,pp.r); 
        pp.w = (p2.w+p4.w)/2; 
        pp.etae = data.ic.calculate_etae_2d(pp.r,pp.psi); 
        pp.We = data.ic.get_We(pp.etae); 
        pp.He = data.ic.get_He(pp.etae); 
        for j = 1:data.plasma.n_electrons
            pp.we{j} = pp.We{j}*pp.r; 
            pp.ne{j} = data.plasma.electrons{j}.n(pp.He{j}-data.plasma.electrons{j}.q*pp.phi);
        end        
        pp.c = data.plasma.cs(pp.ne);
    end
    
    % Slope terms for lines p,m    
    Lpp = dimagno.solver.charline_p(p2.z,p2.r,pp.u,pp.v,pp.c);
    Lp = dimagno.straightline(Lpp);
       
    % Properties at po
    if iter == 0
        po = p3;
    else
        po.u = (p3.u + p4.u)/2; po.v = (p3.v + p4.v)/2;        
    end
      
    %% Obtain position of point 4 and line parameters tp (NEW 20120104)
     
    % searches the intersection between line element Lp and streamline given by psi0
    noborderfound = 1; % Stays =1 if no convergence is reached
    tp = 0; % initial solution (arc length parameter)
    for i=1:data.dimagno.max_iter
        r = Lp.y(tp); % The abs avoids r<0 if iterations oscillate strongly
        [psi,Bz,Br] = data.field.field_2d(Lp.x(tp),r);
        delta = r^2*me_(j_)*We_(j_) - De_(j_) - psi;
        derivative = r*(Br*Lp.A + (2*me_(j_)*We_(j_) - Bz)*Lp.C);
        tp = tp - delta/derivative;
        if abs(delta)+abs(delta/derivative) < data.dimagno.tol
            [Sz,Sr] = data.ic.electron_slope_2d(Bz,Br,num2cell(We_));
            Lw.A = Sz{j_}; % Use only the first species!! 20120411
            Lw.C = Sr{j_};
            Lw.B = Lp.x(tp);
            Lw.D = Lp.y(tp);
            noborderfound = 0;
            break;
        end
    end
    
    if noborderfound == 1
        error([mfilename,':NoBorderFound'],['dMoC border point failed for: z3=',num2str(p3.z),' r3=',num2str(p3.r),' z2=',num2str(p2.z),' r2=',num2str(p2.r),': no border found\n']);        
    end
    
    p4.z = Lw.B; p4.r = Lw.D;
    [p4.psi,p4.Bz,p4.Br] = data.field.field_2d(p4.z,p4.r);
    
    %% canonical momentum equation on charline o (needs nothing more) --> p4.w     
    if me_(j_) == 0
        if data.field.is_axi
            p4.w = p3.w*p3.r/p4.r;
        else
            p4.w = p3.w;
        end
    else
        if data.field.is_axi
            p4.w = (p3.r*p3.w + q_m*(p3.psi - p4.psi)*sim_params.is_magnetic_force_on_ions_on )/p4.r;
        else
            p4.w = p3.w + q_m*(p3.psi - p4.psi)*sim_params.is_magnetic_force_on_ions_on;
        end
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
   
    % System Matrix
    MC = [pp.v, -pp.u, -q_m*sqrt(pp.u^2+pp.v^2-pp.c^2)/pp.c; % char. p
          Lw.C, -Lw.A,                                    0; % border condition
          po.u,  po.v,                                 q_m]; % char. o
    
    % Column vector containing values at previous points
    Gi(1,1) = MC(1,:)*[p2.u;p2.v;p2.phi];
    Gi(2,1) = 0;
    Gi(3,1) = MC(3,:)*[p3.u;p3.v;p3.phi];
    
    % Column vector containing forcing terms (multiplied by the line parameter increments from points 2,1,3, respectively
    Fi(1,1) = ((Lp.A*pp.v-Lp.C*pp.u)*fp + Lp.C*gp - Lp.A*hp) * tp;
    Fi(2,1) = 0;
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
    error([mfilename,':NOCONV'],'Border convergence was not achieved for p4:[z,r] = [%G,%G]',p4.z,p4.r);
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
