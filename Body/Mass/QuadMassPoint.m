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
classdef QuadMassPoint < MassPoint
    % Class used to help compute mass moments of interia
    %
    % A QuadMassPoint is a mass point constructed from a 6 sided
    % trapazoiadal volume. The class computes the volume and centriod using
    % the divergence theorem.
   
     methods
         
        function [mp] = QuadMassPoint(rho, verts)
            % first four points are the 'bottom face' ordered 1,2,3,4 with
            % points in a counter clockwise direction so that the normal
            % points out
            % second four are the 'top face' where the points, 5,6,7,8 are
            % oriented so that 5 is above 1, 6 is above 2 etc.
            
            % uses the divergence theorem to find the volume
            faceInds = [1 2 3 4; 1 4 8 5; 4 3 7 8; 3 2 6 7; 2 1 5 6; 8 7 6 5];
            vol = 0;
            vcent = [0 0 0];
            
            for n = 1:6
                pan = verts(faceInds(n,:), :);
                [norm, ar, cent] = Panel.Properties(pan);
                Fvol = cent(3);
                Fcent = 0.5*cent.^2;
                
                vol = vol + Fvol*norm(3)*ar; % F*n*dA
                vcent = vcent + Fcent.*norm*ar; % Computing 3 divergence 
                                            % theorems for each direction 
                                            % x,y,z, which is why .* is 
                                            % used - represents the dot 
                                            % product in each direction 
            end
            
            vcent = vcent./vol;
            
            mp = mp@MassPoint(rho, vol, vcent);            
        end
     end
end