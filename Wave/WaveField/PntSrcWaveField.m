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
classdef PntSrcWaveField < CirWaveField
    % A wave field provides information (pressure, elevation, velocity,
    % significant wave height, spectrums) at discrete points in space
    % specified on either a grid (array) or at arbitrary (x,y,z) points.  
    %
    % Constructs a wave field of wave radiating from a point source at a
    % specified point on the surface.  The point source can be a heaving,
    % surging or swaying source.

    methods
        % Constructor
        function [wf] = PntSrcWaveField(type, org, rho, waves, isarray, varargin)
            N = waves.Count;
            coefs = cell(N,1);
            
            for n = 1:N
                a = waves.A(n);
                epsilon = waves.Epsilon(n);
                A = a*exp(1i*epsilon);
                switch type
                    case 'Heave'
                        coefs{n} = A;
                    case 'Surge'
                        A1 = 0.5*A;
                        coefs{n} = [-A1, 0, A1];
                    case 'Sway'
                        A1 = -1i*0.5*A;
                        coefs{n} = [A1, 0, A1];
                end
            end               
            
            waves2 = CirWaves('Out', org, coefs, waves.T, waves.H);

            wf = wf@CirWaveField(rho, waves2, isarray, varargin{:});
        end
    end
end