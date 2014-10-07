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
classdef PlaneWaves < IWaves
    
    properties (Dependent)
       IsIncident;      % Indicates whether it is an incident wave.
       IsPlane;         % Indicates whether it is a long-crested plane wave  
       Mlim;            % Truncation values of the circular coefficient for the circular modes. 
       L;               % Truncation values of the circular coefficients for the evanescent modes. 
       KochinFunc;      % The Kochin function is the far-field amplitude as a function of direction.
    end
    
    methods
        
        function [wav] = PlaneWaves(a, t, beta, h, varargin)
            % Constructor
            wav.nT = -1;
            if (nargin ~= 0)
                wav.init(a, t, beta, h, varargin{:});
            end
        end
        
        function [ii] = get.IsIncident(wav)
            % Indicates whether it is an incident wave.
            ii = true;
        end
        
        function [ip] =  get.IsPlane(wav)
            % Indicates whether it is a long-crested plane wave  
            ip = true;
        end
        
        function [mlim] = get.Mlim(wav)
            % Truncation values of the circular coefficient for the circular modes. 
            mlim = Inf*ones(wav.nT, 1);
        end
        
        function [l] = get.L(wav)
            % Truncation values of the circular coefficients for the evanescent modes. 
            l = zeros(wav.nT, 1);
        end
        
        function [koch] = get.KochinFunc(wav)
            % The Kochin function is the far-field amplitude as a function of direction.
            koch = cell(wav.nT, 1); % TODO: get this right
        end
        
        function [as] = IncAmps(wav, M, varargin)
            % Circular coefficients describing the wave.
            
            iwav = 0;
            if (~isempty(varargin))
                pos = varargin{1};
                if (length(varargin) > 1)
                    iwav = varargin{2};
                end
            else
                pos = [0 0];
            end
                
            if (iwav == 0)
                startn = 1;
                stopn = wav.nT;
                as = cell(wav.nT, 1);
            else
                startn = iwav;
                stopn = startn;
            end
            
            for n = startn:stopn
                val = pos(1)*cos(wav.beta(n)) + pos(2)*sin(wav.beta(n));
                exps = zeros(2*M+1,1);
                for m = -M:M
                    exps(m+M+1) = exp(-1i*m*(wav.beta(n) + pi/2));
                end
                A0 = wav.a(n)*exp(1i*wav.epsilon(n))*exp(-1i*wav.k(n)*val);
                
                if (iwav == 0)
                    as{n} = A0*exps;
                else
                    as = A0*exps;
                end
            end
        end
    end
    
end
