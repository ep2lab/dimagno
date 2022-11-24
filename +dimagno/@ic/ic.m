%{
This class provides the methods to query He, De, We and their etae
derivatives, for an initial condition at the z = 0 plane. It also
provides the methods to obtain eta and the electron slope.

Objects of this class need to be initialized during the preprocessing 
with some library vectors: r_,M_,w_,phi_,ne_; and the field and plasma
objects. From them, the library vectors for He, etc, are calculated
internally.  

AFTER CHANGING any of the input vectors, h.recalculate_libraries must be
run.

The class also initializes a front F with method F = initial_front(h,F),
which  is usually called from the preprocessor. The velocities u, v are
directed along B lines. 

The query functions require that r_ be equispaced.
The first element in r_ needs to be on the axis, the latest on the
outermost plasma line. 

MMM20130201
%}
classdef ic < handle & hgsetget
%----------------------------------------------------------------------
    properties % library vectors (given as input to calculate the IC)
        r_; % MUST be equispaced, as produced by linspace!! And ascending from 0 upwards!
        M_; 
        w_;        
        phi_;
        ne_; % a cell array, for when there are more than one electron species
    end 
    properties
        field; % Simply the handles to field and plasma objects
        plasma; 
    end
%----------------------------------------------------------------------
    properties (SetAccess = protected) % Library vectors, calculated with "recalculate_libraries"
        i_; % just 1,2,3... useful array for indexing
        He_;
        De_;
        We_;
        dHe_dr_; % r is eta (the variable that maps the electron streamtubes)
        dDe_dr_;
        dWe_dr_;
        n_r_; % number of points in the library vectors (updated upon updating r_)  
        r_max; % maximum radius at initial condition
    end 
%----------------------------------------------------------------------    
    methods % object preparation actions 
        function h = field_has_changed(h,new_field) 
            % takes the necessary steps when magfield is new (common IC
            % interface)
            h.field = new_field;
            h.recalculate_libraries;
        end
        function h = recalculate_libraries(h) 
            h.He_ = h.calculate_He_;
            h.De_ = h.calculate_De_;
            h.We_ = h.calculate_We_;
            h.dHe_dr_ = h.calculate_dHe_dr_;
            h.dDe_dr_ = h.calculate_dDe_dr_;
            h.dWe_dr_ = h.calculate_dWe_dr_;
            h.i_ = 1:h.n_r_;
        end
    end
%----------------------------------------------------------------------
    methods % Preprocessor actions
        function F = initial_front(h,F) 
            r = F.r_;
            %allocate             
            F.u_ = r.*0;
            F.v_ = r.*0;
            F.w_ = r.*0;
            F.phi_ = r.*0;
            F.c_ = r.*0;            
            F.psi_ = r.*0;
            F.Bz_ = r.*0;
            F.Br_ = r.*0;
            for i = 1: h.plasma.n_electrons
                F.ne_{i} = r.*0;
                F.Te_{i} = r.*0;
                F.we_{i} = r.*0;
                F.etae_{i} = r.*0;
            end
            % calculate this to use the last point slope (wall) below
            [F.psi_(end),F.Bz_(end),F.Br_(end)] = h.field.field_2d(0,r(end));
            distance_to_vertex = F.Bz_(end)/F.Br_(end);            
            % calculate rest of points
            for j = 1:length(r)
                for i = 1: h.plasma.n_electrons
                    temp{i} = h.qinterp1(h.ne_{i},r(j)); % temp
                    F.ne_{i}(j) = temp{i};
                    F.Te_{i}(j) = h.plasma.electrons{i}.T(temp{i});
                    F.we_{i}(j) = h.qinterp1(h.We_{i},r(j))*r(j);
                    F.etae_{i}(j) = r(j);
                end
                F.c_(j) = h.plasma.cs(temp);
                M = h.qinterp1(h.M_,r(j));
                [F.psi_(j),F.Bz_(j),F.Br_(j)] = h.field.field_2d(0,r(j));
                B = sqrt(F.Bz_(j)^2+F.Br_(j)^2);
                F.u_(j) = M*F.c_(j)*cos(atan2(r(j),distance_to_vertex));
                F.v_(j) = M*F.c_(j)*sin(atan2(r(j),distance_to_vertex));
                F.w_(j) = h.qinterp1(h.w_,r(j));
                F.phi_(j) = h.qinterp1(h.phi_,r(j));
            end
        end        
    end
%----------------------------------------------------------------------
    methods % Query libraries 
        function He = get_He(h,etae) % eta is r0 here
            if ~iscell(etae)
                etae = {etae}; % in case of only one species, when only one density array is provided, it is not necessary to give it like a cell
            end
            for i = 1:h.plasma.n_electrons
                He{i} = h.qinterp1(h.He_{i},etae{i});
            end
        end
        function De = get_De(h,etae) % eta is r0 here
            if ~iscell(etae)
                etae = {etae}; % in case of only one species, when only one density array is provided, it is not necessary to give it like a cell
            end
            for i = 1:h.plasma.n_electrons
                De{i} = h.qinterp1(h.De_{i},etae{i});
            end
        end
        function We = get_We(h,etae) % eta is r0 here
            if ~iscell(etae)
                etae = {etae}; % in case of only one species, when only one density array is provided, it is not necessary to give it like a cell
            end
            for i = 1:h.plasma.n_electrons
                We{i} = h.qinterp1(h.We_{i},etae{i});
            end
        end
        function dHe_detae = get_dHe_detae(h,etae) % eta is r0 here
            if ~iscell(etae)
                etae = {etae}; % in case of only one species, when only one density array is provided, it is not necessary to give it like a cell
            end
            for i = 1:h.plasma.n_electrons
                dHe_detae{i} = h.qinterp1(h.dHe_dr_{i},etae{i});
            end
        end
        function dDe_detae = get_dDe_detae(h,etae) % eta is r0 here
            if ~iscell(etae)
                etae = {etae}; % in case of only one species, when only one density array is provided, it is not necessary to give it like a cell
            end
            for i = 1:h.plasma.n_electrons
                dDe_detae{i} = h.qinterp1(h.dDe_dr_{i},etae{i});
            end
        end
        function dWe_detae = get_dWe_detae(h,etae) % eta is r0 here
            if ~iscell(etae)
                etae = {etae}; % in case of only one species, when only one density array is provided, it is not necessary to give it like a cell
            end
            for i = 1:h.plasma.n_electrons
                dWe_detae{i} = h.qinterp1(h.dWe_dr_{i},etae{i});
            end
        end
    end
    methods % Calculate eta, slopes and grad_He
        function etae = calculate_etae_2d(h,r,psi)
            for i = 1:h.plasma.n_electrons
                if h.plasma.electrons{i}.m == 0
                    library = h.De_{i}/h.plasma.electrons{i}.q;
                    if psi> max(library) % To extrapolate with constant values if need be.
                        etae{i} = h.r_max;
                        continue; % cycle
                    elseif psi < min(library)                        
                        etae{i} = 0;
                        continue;
                    end
                    etae{i} = interp1(library,h.r_,psi);
                else
                    % This method only works if (De - me*r^2*We) is a monotonic function
                    % for the given r. So it eventually breaks downstream! 
                    % MMM: the extrapolation needs to be implemented in
                    % this case too
                    for j = 1:numel(psi)
                        if h.field.is_axi
                            etae{i}(j) = interp1((h.De_{i}-r(j)^2*h.plasma.electrons{i}.m.*h.We_{i})/h.plasma.electrons{i}.q,h.r_,psi(j));
                        else
                            etae{i}(j) = interp1((h.De_{i}-r(j)*h.plasma.electrons{i}.m.*h.We_{i})/h.plasma.electrons{i}.q,h.r_,psi(j));                     % planar case
                        end
                    end
                end                 
            end                  
        end
        function [Sez,Ser] = electron_slope_2d(h,Bz,Br,We)
            % Returns Sez,Ser: arrays with one element for each electron
            % species (length of cellarray We), which put together give ure/uze
            % = Ser/Sez. Outputs are cell arrays, too
            Ser = We; Sez = We; % allocate
            for i = 1:h.plasma.n_electrons
                Ser{i} = Br;
                Sez{i} = Bz-2*h.plasma.electrons{i}.m*We{i};
            end
        end    
        function [dHe_dz, dHe_dr] = calculate_grad_He_2d(h,r,psi,Bz,Br)
            % It takes the field components already, so no need to recalculate
            % them.              
            [dHe_deta,dDe_deta,dWe_deta,We] = h.calculate_dHe_deta_dDe_deta_dWe_deta_We_2d(r,psi);
            for i = 1:h.plasma.n_electrons
                m = h.plasma.electrons{i}.m;
                if h.field.is_axi == 0
                    m = 0; % in planar cases, remove this contribution
                end
                q = h.plasma.electrons{i}.q;
                part = dHe_deta{i}/(dDe_deta{i}-m*r^2*dWe_deta{i});
                dHe_dz{i} = -q*r*Br * part;
                dHe_dr{i} = (q*r*Bz+2*r*m*We{i}) * part;
            end
        end
    end 
%----------------------------------------------------------------------
    methods % BASIC
        function set.r_(h,v) % simply sets n_r_ at the same time
            h.r_ = v;
            h.n_r_ = length(h.r_);
            h.r_max = h.r_(h.n_r_); % maximum value of r
        end
        function set.ne_(h,v)
            if ~iscell(v) % transform it to cell if not cell
                v = {v};
            end
            h.ne_ = v(:); % force column
        end
        function h = ic(varargin)
            % Parse input
            p = inputParser; 
            p.addParameter('field',magnetic_field.loop_2d,@(x)isa(x,'magnetic_field.element_2d')); % Field object
            p.addParameter('plasma',fluid_plasma.plasma,@(x)isa(x,'fluid_plasma.plasma')); % Plasma object

            % Find out how many species does plasma have
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            n_electrons = p.Results.plasma.n_electrons;
            p.KeepUnmatched = false;

            % Continue the validation     
            p.addParameter('r_',linspace(0,1,1000),@isnumeric);
            p.addParameter('M_',linspace(0,1,1000).*0+1.01,@isnumeric);
            p.addParameter('w_',linspace(0,1,1000).*0,@isnumeric);
            p.addParameter('phi_',linspace(0,1,1000).*0,@isnumeric);
                default_ne = exp(-(linspace(0,1,1000).^2)*3*log(10)); 
                default_necell(1:n_electrons) = {default_ne};
            p.addParameter('ne_',default_necell,@iscell); 
            
            p.parse(varargin{:}); % check all, and assign defaults to p.Results as needed.
            fn = fieldnames(p.Results);
            for i = 1:length(fn)
                h.(fn{i}) = p.Results.(fn{i});
            end            
            % Calculate libraries
            h.recalculate_libraries;
        end           
    end
    methods (Hidden = true, Access = protected)
        function Yi = qinterp1(h,Y,eta) % ONLY USABLE BECAUSE r_ is linspace-like! Y is the library to use
            % eta must be scalar!
            t = h.r_max; % maximum value of r
            % Calculate eta in index space.
            s = eta*(h.n_r_-1)/t + 1;  
            % Linear interpolation method, extrapolating last value outside
            s = min(h.n_r_,max(1,s));
            s0 = floor(s);
            s1 = ceil(s);
            ds = s - s0; % excess or increment from s0
            Yi = Y(s0).*(1-ds) + Y(s1).*ds;
        end        
        function [dHe_deta,dDe_deta,dWe_deta,We] = calculate_dHe_deta_dDe_deta_dWe_deta_We_2d(h,r,psi)
            % This function tries to save time by calculating all at the same
            % time. It uses a variant of qinterp1, implicitly
            % First, calculate in index space the position to take from the
            % libraries
            for j = 1:h.plasma.n_electrons
                if h.plasma.electrons{j}.m == 0
                    i__  = interp1(h.De_{j}/h.plasma.electrons{j}.q,h.i_,psi,'linear','extrap');
                else
                    % This method only works if (De - me*r^2*We) is a monotonic function
                    % for the given r. So it eventually breaks downstream!
                    % PSI MUST BE SCALAR!! 
                    if h.field.is_axi
                        i__  = interp1((h.De_{j}-r^2*h.plasma.electrons{j}.m.*h.We_{j})/h.plasma.electrons{j}.q,h.i_,psi,'linear','extrap');
                    else
                        i__  = interp1((h.De_{j}-r*h.plasma.electrons{j}.m.*h.We_{j})/h.plasma.electrons{j}.q,h.i_,psi,'linear','extrap');    % planar case
                    end                    
                end                
                % Linear interpolation method
                si_ = floor(i__); % s is now the floor of eta % this is just an optimization
                ti_ = max(1,min(h.n_r_,si_+1));  % limit them from 1...n_r_. This serves to EXTRAPOLATE directly
                si_ = max(1,min(h.n_r_,si_));  % limit them from 1...n_r_. This serves to EXTRAPOLATE directly
                i__ = i__-si_;                 % increment (saved on eta to save space)
                dHe_deta{j} = h.dHe_dr_{j}(si_).*(1-i__) + h.dHe_dr_{j}(ti_).*i__;
                dDe_deta{j} = h.dDe_dr_{j}(si_).*(1-i__) + h.dDe_dr_{j}(ti_).*i__;
                dWe_deta{j} = h.dWe_dr_{j}(si_).*(1-i__) + h.dWe_dr_{j}(ti_).*i__;
                We{j} = h.We_{j}(si_).*(1-i__) + h.We_{j}(ti_).*i__;
            end                  
        end
        function we_ = calculate_we_(h)
            [~,Bz_,~] = h.field.field_2d(h.r_.*0,h.r_);
            dhe_dr_ = h.calculate_dhe_dr_;
            dphi_dr_ = h.calculate_dphi_dr_;
            for i = 1:h.plasma.n_electrons
                e = h.plasma.electrons{i};
                if e.m ~= 0 && h.field.is_axi % if not axi, inertia plays no role
                    we_{i} = -e.q*Bz_.* h.r_ ./ (2.*e.m);
                    we_{i} = we_{i} - sqrt(e.q^2*Bz_.^2.* h.r_.^2 ./ (4.*e.m^2) + (dhe_dr_{i} + e.q*dphi_dr_).* h.r_ / e.m);
                    
                else
                    we_{i} =  ( dhe_dr_{i} + e.q *dphi_dr_) ./ (e.q * Bz_); % when no inertia
                end
            end
        end
        function He_ = calculate_He_(h)
            we = h.calculate_we_;
            for i = 1:h.plasma.n_electrons
                e = h.plasma.electrons{i};
                He_{i} = e.h(h.ne_{i}) + e.q*h.phi_ + e.m/2*we{i}.^2; % we can keep this contribution to He even if planar
            end
        end
        function De_ = calculate_De_(h)
            we = h.calculate_we_;
            for i = 1:h.plasma.n_electrons
                e = h.plasma.electrons{i};
                if  h.field.is_axi 
                    De_{i} = h.r_*e.m.*we{i} + e.q*h.calculate_psi_;
                else
                    De_{i} = e.m.*we{i} + e.q*h.calculate_psi_;
                end
            end
        end
        function We_ = calculate_We_(h) % electron azimuthal frequency
            we = h.calculate_we_;
            [~,Bz,~] = h.field.field_2d(h.r_(1).*0,h.r_(1));
            for i = 1:h.plasma.n_electrons
                e = h.plasma.electrons{i};
                We_{i} = we{i}./h.r_;
                % Value at the origin
                he1 = e.h(h.ne_{i}(1));
                he2 = e.h(h.ne_{i}(2));
                We_{i}(1) = 1./(e.q*Bz) * 2*(he2-he1+h.phi_(2)-h.phi_(1))/h.r_(2)^2;
            end            
        end
        function dHe_dr_ = calculate_dHe_dr_(h)
            for i = 1:h.plasma.n_electrons
                dHe_dr_{i} = utilities.gradient(h.r_.',h.He_{i}.').';  
            end
        end
        function dDe_dr_ = calculate_dDe_dr_(h)
            for i = 1:h.plasma.n_electrons
                dDe_dr_{i} = utilities.gradient(h.r_.',h.De_{i}.').'; 
            end
        end
        function dWe_dr_ = calculate_dWe_dr_(h)
            for i = 1:h.plasma.n_electrons
                dWe_dr_{i} = utilities.gradient(h.r_.',h.We_{i}.').';
            end
        end
        function dhe_dr_ = calculate_dhe_dr_(h)
            for i = 1:h.plasma.n_electrons
                e = h.plasma.electrons{i};
                he_ = e.h(h.ne_{i});
                dhe_dr_{i} = utilities.gradient(h.r_.',he_.').';
            end
        end
        function dphi_dr_ = calculate_dphi_dr_(h)
            dphi_dr_ = utilities.gradient(h.r_.',h.phi_.').'; 
        end
        function psi_ = calculate_psi_(h)
            psi_ = h.field.field_2d(h.r_*0,h.r_);            
        end
    end
%----------------------------------------------------------------------
end
