%{ 
mwave - A water wave and wave energy converter computation package 
Copyright (C) 2014  Cameron McNatt

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contributors:
    C. McNatt
%}
classdef MassPoint
    % Class used to help compute mass moments of interia
    
    properties (Access = protected)
        rho;
        vol;
        pos;
    end

    properties (Dependent)
        Rho;
        Mass;
        Vol; 
        Pos;
    end
    
     methods
         
        function [mp] = MassPoint(varargin)
            %Constructor
            if (isempty(varargin))
                mp.rho = 0;
                mp.vol = 0;
                mp.pos = [0 0 0];
            elseif (length(varargin) == 1)
                copymp = varargin{1};
                if (~isa(copymp, 'MassPoint'));
                    error('To copy a mass point, input must be a MassPoint')
                end
                mp.rho = copymp.rho;
                mp.vol = copymp.vol;
                mp.pos = copymp.pos;
            elseif (length(varargin) == 3)
                rh = varargin{1};
                vo = varargin{2};
                ps = varargin{3};
                if (~isscalar(rh))
                    error('rho input must be scalar');
                end

                if (~isscalar(vo))
                    error('vol input must be scalar');
                end

                if (rh <= 0)
                    error('rho input must be positive');
                end

                if (vo <= 0)
                    error('vol input must be positive');
                end

                [row, col] = size(ps);
                if((row ~= 1) || (col ~= 3))
                    error('pos must be a 1x3 vector');
                end

                mp.rho = rh;
                mp.vol = vo;
                mp.pos = ps;
            else
                error('MassPoint constructor argument not recognized');
            end
        end
        
        function [rh] = get.Rho(mp)
            % The mass density
            rh = mp.rho;
        end
        function [mp] = set.Rho(mp, rh)
            % The mass density
            if (~isscalar(rh))
                error('Rho input must be scalar');
            end
                       
            if (rh <= 0)
                error('Rho input must be positive');
            end
            
            mp.rho = rh;
        end
        
        function [vo] = get.Vol(mp)
            % The volume of mass point
            vo = mp.vol;
        end
        function [mp] = set.Vol(mp, vo)
            % The volume
            if (~isscalar(vo))
                error('Vol input must be scalar');
            end
                       
            if (vo <= 0)
                error('Vol input must be positive');
            end
            
            mp.vol = vo;
        end
        
        function [m] = get.Mass(mp)
            % The mass of the poit
            m = mp.rho*mp.vol;
        end
        
        function [p] = get.Pos(mp)
            % The point position
            p = mp.pos;
        end
        function [mp] = set.Pos(mp, p)
            % The point position
            [row, col] = size(p);
            if((row ~= 1) || (col ~= 3))
                error('pos must be a 1x3 vector');
            end
            mp.pos = p;
        end
     end
end