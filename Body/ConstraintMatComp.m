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
classdef ConstraintMatComp
    
    methods (Static)
        function [P] = HingedBodies(bods, hins, varargin)
            % Inputs:
            %   bods = N x 3 matrix of {x, y, z} body coordinates (location
            %       of the body in global coordinates), where N is the 
            %       number of bodies
            %   hins = (N-1) x 3 matrix of {x, y, z} hinge coordinates. The
            %       y-coordinate is not really necessary as it's a hinge 
            %       about an axis parallel to the y-axis
            %
            % Optional Inputs:
            %   'Origin', org = org is a 1 x 3 vector indicating the origin
            %       of the composite body. If the optional 'Origin'
            %       argument is not provided, the default is is the body
            %       coordinates of body 1
            %   'Planar': compute the constraint matrix for 3 DOF x 3 DOF
            %   planar motions, where the DOF are surge, heave, and pitch
            %
            % Returns:
            %   P = the velocity transformation matrix 
            %       (not PT, i.e. the transpose) 
            
            [opts, args] = checkOptions({{'Origin', 1}, {'Planar'}}, varargin);
            
            if (opts(1))
                org = args{1};
            else
                org = bods(1,:);
            end
            
            planar = opts(2);
            
            Nbod = size(bods, 1);
            Nhin = size(hins, 1);
            
            if (Nhin ~= (Nbod - 1))
                error(['The number of hinges must be one less than the '...
                'number of bodies']);
            end
            
            if (Nbod > 1)
                if ((size(bods, 2) ~= 3) || (size(hins, 2) ~= 3))
                    error(['The body coordinates and hinge coordinates '...
                        'must have x,y,z locations']);
                end
            end
            
            PT = zeros(6*Nbod, 6 + Nhin);
            
            sR = zeros(Nbod, 3);
            sL = zeros(Nbod, 3);
            
            for n = 1:Nbod
                if (n < Nbod)
                    sR(n,:) = hins(n,:) - bods(n,:);
                end
                
                if (n > 1)
                    sL(n,:) = hins(n-1,:) - bods(n,:);
                end
            end
            
            sO = bods(1,:) - org;
              
            for n = 1:Nbod
                PTn = zeros(6, 6 + Nhin);
                
                % Identity matrices
                PTn(1:3, 1:3) = eye(3);
                PTn(4:6, 4:6) = eye(3);
                                
                % cross-product matrix
                svect = -sO;
                for m = 2:n
                    svect = svect - sR(m-1,:) + sL(m,:);
                end
                Sx = skewMat(svect);
                PTn(1:3,4:6) = Sx;
                
                for o = 2:n
                    svect = [-sL(n,3), 0, sL(n,1)];
                
                    for m = n:-1:(o+1)
                        svect = svect + [-sL(m-1,3), 0, sL(m-1,1)] ...
                            + [sR(m-1,3), 0, -sR(m-1,1)];
                    end
                    
                    PTn(1:3,o+5) = svect';
                end

                % row of 1's in pitch for flex modes
                PTn(5,7:(5+n)) = ones(1,n-1);
                                
                istart = (n - 1)*6 + 1;
                PT(istart:(istart + 5), :) = PTn;
            end
            
            P = PT.';
            
            if planar
%                 [M, N] = size(P);
%                 P = P(1:2:M, 1:2:N);
                P(2:2:6,:) =[]; % remove constrained out of plane
                P(:,2:2:end) =[]; % remove body out of plane
            end
        end
        
        function [P] = FixedBodies(bods, varargin)
            Nb = size(bods,1);
            Nhin = Nb-1;
            P0 = ConstraintMatComp.HingedBodies(bods, zeros(Nhin,3), varargin{:});

            P = P0(1:(end-Nhin),:);
        end
        
        function [P] = RelativeHeave(bods, varargin)
            % Inputs:
            %   bods = 2 x 3 matrix of {x, y, z} body coordinates (location
            %       of the body in global coordinates), where N = 2 is the 
            %       number of bodies
            %
            % Optional Inputs:
            %   'Origin', org = org is a 1 x 3 vector indicating the origin
            %       of the composite body. If the optional 'Origin'
            %       argument is not provided, the default is is the body
            %       coordinates of body 1
            %   'Planar': compute the constraint matrix for 3 DOF x 3 DOF
            %   planar motions, where the DOF are surge, heave, and pitch
            %
            % Returns:
            %   P = the velocity transformation matrix 
            %       (not PT, i.e. the transpose) 
            
            [opts, args] = checkOptions({{'Origin', 1}, {'Planar'}}, varargin);
            
            if (opts(1))
                org = args{1};
            else
                org = bods(1,:);
            end
            
            planar = opts(2);
            
            % sO is vector from origin to body 1
            sO = bods(1,:) - org;
            
            % s12 is vector from body 1 to body 2
            s12 = bods(2,:) - bods(1,:);
                        
            % body 1 has same angles as origin
            % body 1 CG is origin + rotation of relative body vector (sO)
            PT1 = zeros(6, 7);
            PT1(1:6, 1:6) = eye(6);
            PT1(1:3, 4:6) = -skewMat(sO);
                                    
            % body 2 has same angles as origin
            % body 2 CG =
            %   origin 
            %   + rotation of relative body vector (sO + s12)
            %   + relative heave vector [0, 0, zh]'
            PT2 = zeros(6, 7);
            PT2(1:6, 1:6) = eye(6);
            PT2(1:3, 4:6) = -skewMat(sO + s12);
            PT2(3, 7) = 1;
            
            % Note that the vectors in the skewMat (sO, and sO + s12) have
            % a negative sign because: R*s = s - skewMat(s)*alpha, where 
            % R is linear the rotation matrix
            % s on the rhs is not part of P bc it is fixed, and alpha is
            % the vector of rotation angles
            
            PT = [PT1; PT2];
            
            P = PT.';
            
            if planar
                [M, N] = size(P);
                P = P(1:2:M, 1:2:N);
            end
        end
    end
end