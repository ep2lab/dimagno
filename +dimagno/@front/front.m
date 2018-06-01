%----------------------------------------------------------------------
%{
This abstract class defines the basic interface for simulation fronts. 
It keeps enough data from the last
calculated points as to be able to continue the integration from there.
 
front should be a surface/line out of the Mach cone of each of its points.

MMM20130221
%}
%----------------------------------------------------------------------
classdef front < handle
    properties % Basic Properties (arrays) in the front (from axis (1) to wall (end))        
        z_
        r_
        u_
        v_
        w_
        phi_
        ne_ % Electron property: cell array of arrays
        c_        
    end
    properties % derived variables (unimportant, but kept)
        psi_
        Bz_
        Br_ 
        Te_ % Electron property: cell array of arrays
        we_ % Electron property: cell array of arrays       
        etae_ % Electron property: cell array of arrays
    end
    properties (Dependent = true)
        M_ % calculated as sqrt(u_.^2+v_.^2)./c_
    end
%----------------------------------------------------------------------
    properties 
        n_points; % number of points
        n_electrons; % number of electron species
    end
%----------------------------------------------------------------------    
    properties (Constant = true)
        basic_variable_names =  {'z','r','u','v','w','phi','c'};
        basic_electron_variable_names =  {'ne'};
        derived_variable_names =  {'psi','Bz','Br'};
        derived_electron_variable_names =  {'Te','we','etae'};
        dependent_variable_names = {'M'};
        basic_variable_names_ =  {'z_','r_','u_','v_','w_','phi_','c_'};
        basic_electron_variable_names_ =  {'ne_'};
        derived_variable_names_ =  {'psi_','Bz_','Br_'};
        derived_electron_variable_names_ =  {'Te_','we_','etae_'};
        dependent_variable_names_ = {'M_'};
    end
%----------------------------------------------------------------------
    methods
        function p = point(h,varargin)
            % Returns a structure with the selected point. Varargin is just the
            % indicing to that point.
            for var = h.basic_variable_names
                p.(var{1}) = h.([var{1},'_'])(varargin{:});                
            end
            for var = h.basic_electron_variable_names
                for j = 1:h.n_electrons
                    p.(var{1}){j} = h.([var{1},'_']){j}(varargin{:});
                end
            end
            for var = h.derived_variable_names
                try
                    p.(var{1}) = h.([var{1},'_'])(varargin{:});                
                end
            end
            for var = h.derived_electron_variable_names
                try
                    for j = 1:h.n_electrons
                        p.(var{1}){j} = h.([var{1},'_']){j}(varargin{:});
                    end
                end
            end
        end
        function h = set.z_(h,v)
            h.z_ = v; 
            h.n_points = numel(h.z_); % we do this on setting "points"
        end
        function h = set.ne_(h,v)
            if ~iscell(v)
                v = {v};
            end
            h.ne_ = v;             
            h.n_electrons = numel(v);
        end
        function h = set.Te_(h,v)
            if ~iscell(v)
                v = {v};
            end
            h.Te_ = v;             
        end
        function h = set.we_(h,v)
            if ~iscell(v)
                v = {v};
            end
            h.we_ = v;             
        end
        function h = set.etae_(h,v)
            if ~iscell(v)
                v = {v};
            end
            h.etae_ = v;             
        end
        function v = get.M_(h)
            v = sqrt(h.u_.^2+h.v_.^2)./h.c_;
        end
    end    
    methods % Interpolators for z r u v w phi ne c (the basic variables ONLY)
        function z  = z(h,sf)
            z = h.qinterp1(h.z_,sf);            
        end
        function r = r(h,sf)
            r = h.qinterp1(h.r_,sf);            
        end
        function u  = u(h,sf)
            u = h.qinterp1(h.u_,sf);            
        end
        function v = v(h,sf)
            v = h.qinterp1(h.v_,sf);            
        end
        function w = w(h,sf)
            w = h.qinterp1(h.w_,sf);            
        end
        function phi = phi(h,sf)
            phi = h.qinterp1(h.phi_,sf);
        end
        function ne = ne(h,sf)
            for i=1:h.n_electrons
                ne{i} = h.qinterp1(h.ne_{i},sf);            
            end
        end
        function c = c(h,sf)
            c = h.qinterp1(h.c_,sf);            
        end
        function [z,r,u,v,w,phi,ne,c] = zruvwphinec(h,sf) % combined interpolator for all vars at a time
            % allocate
            z = NaN*sf; r = z; u = z; v = z; w = z; phi = z; c = z;
            % Throw away already those sf not between 0 and 1 (and not Inf nor NaN)
            valid = (sf>=0 & sf<=1 & isfinite(sf));
            sf = sf(valid);
            % Calculate sf in index space.
            sf = sf*(h.n_points-1)+1;            
            % Linear interpolation method
            fsf = floor(sf);          % floor
            csf = ceil(sf);           % ceiling (needs both for when it is on the border)
            sf = sf-fsf;              % increment (saved on xi to save space)
            z(valid) = h.z_(fsf).*(1-sf) + h.z_(csf).*sf;
            r(valid) = h.r_(fsf).*(1-sf) + h.r_(csf).*sf;
            u(valid) = h.u_(fsf).*(1-sf) + h.u_(csf).*sf;
            v(valid) = h.v_(fsf).*(1-sf) + h.v_(csf).*sf;
            w(valid) = h.w_(fsf).*(1-sf) + h.w_(csf).*sf;
            phi(valid) = h.phi_(fsf).*(1-sf) + h.phi_(csf).*sf;
            for i=1:h.n_electrons
                ne{i}(valid) = h.ne_{i}(fsf).*(1-sf) + h.ne_{i}(csf).*sf;         
            end
            c(valid) = h.c_(fsf).*(1-sf) + h.c_(csf).*sf;
        end
    end
    methods 
        function lh = plot_front(h,varargin)
            lh = line(h.z_,h.r_,...
                 'Color','k','Marker','.',...
                 'MarkerFaceColor','r','MarkerEdgeColor','r',...
                 'DisplayName',['Points of front: ',inputname(1)],varargin{:});
        end
        function lh = plot_front_var(h,var,varargin) 
            % When plotting electron properties (cells), you can pass an integer
            % to select which species as the first varargin.
            thing = h.([var,'_']);
            if iscell(thing)
                if ~isempty(varargin) && isnumeric(varargin{1})
                    ind = varargin{1}; % electron species to plot
                    varargin = varargin(2:end);
                else 
                    ind = 1; % just plot first electron species
                end
                thing = thing{ind}; 
                var = [var,'{',num2str(ind),'}'];
            end            
            lh = line(h.z_,h.r_,thing,...
                'Color','k','Marker','.',...
                'MarkerFaceColor','r','MarkerEdgeColor','r',...
                'DisplayName',[var,' of front: ',inputname(1)],varargin{:});
        end        
        function h = setpoint(h,p,pos)
            % this method allows to give a structure and take from its fields
            % the property values and store them in the front arrays in
            % position pos. Pos can also be prepend and append.           
            if strcmp(pos,'prepend') && ~isempty(h.z_) % this tells setpoint to displace all arrays one position
                for var = h.basic_variable_names
                    h.([var{1},'_'])(2:end+1,1) = h.([var{1},'_']);
                end
                for var = h.derived_variable_names
                    if ~isempty(h.([var{1},'_']))
                        h.([var{1},'_'])(2:end+1,1) = h.([var{1},'_']);
                    end
                end                
                for var = h.basic_electron_variable_names
                    for i=1:h.n_electrons
                        h.([var{1},'_']){i}(2:end+1,1) = h.([var{1},'_']){i};
                    end
                end
                for var = h.derived_electron_variable_names
                    if ~isempty(h.([var{1},'_']))
                        for i=1:h.n_electrons
                            h.([var{1},'_']){i}(2:end+1,1) = h.([var{1},'_']){i};
                        end
                    end
                end      
                pos = 1;
            elseif strcmp(pos,'append')
                pos = h.n_points +1;
            end          
            % Basic Properties
            for var = h.basic_variable_names 
                h.([var{1},'_'])(pos,1) = p.(var{1});
            end
            for var = h.basic_electron_variable_names
                for i=1:length(p.ne)
                    h.([var{1},'_']){i}(pos,1) = p.(var{1}){i};
                end
            end
            % Derived properties
            for var = h.derived_variable_names
                if isfield(p,var{1})
                    h.([var{1},'_'])(pos,1) = p.(var{1});
                end
            end
            for var = h.derived_electron_variable_names
                if isfield(p,var{1})
                    for i=1:length(p.ne)
                        h.([var{1},'_']){i}(pos,1) = p.(var{1}){i};
                    end
                end
            end
        end
        function h = double_front(h)
            % This method interpolates within the front, to produce more
            % points in it.            
            % Space all points with one space
            nmax = h.n_points;
            for i=nmax:-1:1 % up to down, not to overwrite
                h.setpoint(h.point(i),2*i-1);
            end            
            % Create and interpolate for all intermediate points
            for i=2:nmax*2-2
                p1 = h.point(i-1);
                p2 = h.point(i+1); 
                for var = h.basic_variable_names
                    p.(var{1}) = (p1.(var{1}) + p2.(var{1}))/2;
                end
                for var = h.derived_variable_names                    
                    try
                        p.(var{1}) = (p1.(var{1}) + p2.(var{1}))/2;
                    end
                end
                for var = h.basic_electron_variable_names
                    for j=1:h.n_electrons
                        p.(var{1}){j} = (p1.(var{1}){j} + p2.(var{1}){j})/2;
                    end
                end
                for var = h.derived_electron_variable_names
                    try
                        for j=1:h.n_electrons
                            p.(var{1}){j} = (p1.(var{1}){j} + p2.(var{1}){j})/2;
                        end
                    end
                end
                h.setpoint(p,i);
            end                
        end
    end
%----------------------------------------------------------------------
    methods (Hidden = true, Access = protected)
        function Yi = qinterp1(h,Y,sf) % ONLY USABLE HERE
            % allocate
            Yi = NaN*sf;
            % Throw away already those sf not between 0 and 1 (and not Inf nor NaN)
            valid = (sf>=0 & sf<=1 & isfinite(sf));
            sf = sf(valid);
            % Calculate sf in index space.
            sf = sf*(h.n_points-1)+1;            
            % Linear interpolation method
            fsf = floor(sf);          % floor
            csf = ceil(sf);           % ceiling (needs both for when it is on the border)
            sf = sf-fsf;              % increment (saved on xi to save space)
            Yi(valid) = Y(fsf).*(1-sf) + Y(csf).*sf;
        end
    end  
%----------------------------------------------------------------------    
    methods
        function h = front(varargin) 
            % Parse input
            p = inputParser; 
            p.addParameter('z_',0); 
            p.addParameter('r_',0); 
            p.addParameter('u_',0);
            p.addParameter('v_',0);
            p.addParameter('w_',0);
            p.addParameter('phi_',0);
            p.addParameter('ne_',{0});
            p.addParameter('c_',0);
            p.addParameter('psi_',0);
            p.addParameter('Bz_',0);
            p.addParameter('Br_',0);
            p.addParameter('Te_',{0});
            p.addParameter('we_',{0});
            p.addParameter('etae_',{0});
 
            p.parse(varargin{:}); % check all, and assign defaults to p.Results as needed.
            fn = fieldnames(p.Results);
            for i = 1:length(fn)
                h.(fn{i}) = p.Results.(fn{i});
            end            
        end
    end
%----------------------------------------------------------------------    
end
