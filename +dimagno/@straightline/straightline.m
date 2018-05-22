%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% #DESCRIPTION: Line structure. Just the structure plus the intersection
% function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% #NOTES: line structures have A,B,C,D so that x = A*t + B; y = C*t + D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% #VERSION HISTORY:
% MMM20110523: First version
% MMM20111130: Add the x and y methods. Add the z dimension, although
% intersection is not yet supported
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% #AUTHOR: Mario Merino Martinex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

classdef straightline < handle
    properties
        A = 1; % Line parameters
        B = 0;
        C = 0;
        D = 0;
        E = 0;
        F = 0;
    end
    methods
        function x = x(h,t)
            % Evaluates x
            if ~exist('t','var')
                t = 0;
            end
            x = h.A.*t + h.B;
        end
        function y = y(h,t)
            % Evaluates  y
            if ~exist('t','var')
                t = 0;
            end
            y = h.C.*t + h.D;
        end
        function z = z(h,t)
            % Evaluates  y
            if ~exist('t','var')
                t = 0;
            end
            z = h.E.*t + h.F;
        end
        function h = normalize_slope(h)
            sq = sqrt(h.A^2+h.C^2+h.E^2);
            h.A = h.A/sq;
            h.C = h.C/sq;
            h.E = h.E/sq;
        end
        function [x,y,to,ta] = intersectionxy(h,a)
            % Only for planar lines (z will be ignored): Calculates the intersection point with another line object
            denom = a.A.*h.C - a.C.*h.A;
            if denom == 0 % lines are parallel
                x = Inf;
                y = Inf;
                to = Inf;
                ta = Inf;
                return;
            end
            to = (a.C.*(h.B-a.B) - a.A.*(h.D-a.D)) ./ denom;
            ta = - (h.C.*(a.B-h.B) - h.A.*(a.D-h.D)) ./ denom;            
            x = h.x(to);
            y = h.y(to);
        end 
        function [d,to] = distancexy(h,x,y)
            % distance (with sign) to the line element of point (x,y). 20121113
            % to is the parameter to the foot of the perpendicular
            x = x-h.B;
            y = y-h.D;
            d = h.A*y - h.C*x; % vector product (it gives the sign: + if to the left of line's vector)
            n = sqrt(h.A^2+h.C^2);
            d = d/n;
            to = x*h.A + y*h.C; % dot product gives the parameter to
            to = to/n^2;
        end 
        function lh = plotline(h,t,varargin)
            % Fast plot of the line (t=0..10)
            if ~exist('t','var')
                t = 10;
            end
            lh =line([h.B,h.B+t*h.A],[h.D,h.D+t*h.C],[h.F,h.F+t*h.E],varargin{:});
        end        
    end
    methods
        function h = straightline(varargin)
            % Parse input
            p = inputParser;
            p.addParameter('A',1); 
            p.addParameter('B',0); 
            p.addParameter('C',0); 
            p.addParameter('D',0); 
            p.addParameter('E',0); 
            p.addParameter('F',0);  
            p.parse(varargin{:}); % check all, and assign defaults to p.Results as needed.
            fn = fieldnames(p.Results);
            for i = 1:length(fn)
                h.(fn{i}) = p.Results.(fn{i});
            end 
        end
    end
    methods (Static = true)
        function h = straightline_from_struct(L)
            % provides a straightline element from a struct with the ABCD fields
            h = straightline;
            h.A = L.A;
            h.B = L.B;
            h.C = L.C;
            h.D = L.D;
            if isfield(L,'E')
                h.E = L.E;
                h.F = L.F;
            end
        end
    end
end
