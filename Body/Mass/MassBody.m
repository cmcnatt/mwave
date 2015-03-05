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
classdef MassBody < handle
    % Class used to help compute mass matrices
    
    properties (Access = private)
        constRho;
        mps;
        modes;
        iscomp
        M;
        mas;
        cg;
    end

    properties (Dependent)
        Rho;
        Modes;
        MassMatrix;
        Mass;
        Cg;
    end
    
     methods
         
        function [mass] = MassBody(massPoints)
            N = length(massPoints);
            m(N,1) = MassPoint;
            for n = 1:N
                if (~isa(massPoints(n), 'MassPoint'));
                    error('All inputs must be MassPoints');
                end
                
                if (n == 1)
                    rho0 = massPoints(1).Rho;
                    mass.constRho = true;
                else
                    if (massPoints(n).Rho ~= rho0)
                        mass.constRho = false;
                    end
                end
                
                m(n) = MassPoint(massPoints(n));
            end
            mass.mps = m;
            mass.modes = ModesOfMotion;
            mass.iscomp = false;
        end
        
        function [rho] = get.Rho(mass)
            % The mass density
            if (mass.constRho)
                rho = mass.mps(1).Rho;
            else
                rho = 'Not constant density';
            end
        end
        function [mass] = set.Rho(mass, rho)
            % The mass density
            if (mass.constRho)
                if (~isscalar(rho))
                    error('Rho input must be scalar');
                end
                       
                if (rho <= 0)
                    error('Rho input must be positive');
                end
                
                rhoOld = mass.mps(1).Rho;
            
                N = length(mass.mps);
                for n = 1:N
                    mass.mps(n).Rho = rho;
                end
                
                if (mass.iscomp)
                    mass.M = rho/rhoOld*mass.M;
                end
            else
                error('Not constant density');
            end
        end
        
        function [mo] = get.Modes(mass)
            % Get the modes of motion
            mo = mass.modes;
        end
        function [mass] = set.Modes(mass, mo)
            % Set the modes of motion
            if (isa(mo, 'ModesOfMotion'))
                mass.modes = mo;
            else
                error('The Modes must be of type ModesOfMotion');
            end
            
            mass.iscomp = false;
        end
        
        function [M_] = get.MassMatrix(mass)
            % The mass matrix
            if (isempty(mass.modes))
                error('Cannot compute mass matrix without modes of motion');
            end
            
            if (~mass.iscomp)
                mass.computeM();
            end
            M_ = mass.M;
        end
        
        function [ma] = get.Mass(mass)
            % The mass
            if (~mass.iscomp)
                mass.computeM();
            end
            ma = mass.mas;
        end
        
        function [c] = get.Cg(mass)
            % The cener of gravity
            if (~mass.iscomp)
                mass.computeM();
            end
            c = mass.cg;
        end
        
        function [] = Translate(mass, vector)
            % Translate the mass
            [row, col] = size(vector);
            if((row ~= 1) || (col ~= 3))
                error('the translation vector must be a 1x3 vector');
            end
            N = length(mass.mps);
            for n = 1:N
                mass.mps(n).Pos = mass.mps(n).Pos + vector;
            end
        end
        
        function [mout] = plus(ma, mb)
            % Overloaded plus operator
            if (~isa(ma, 'Mass') || ~isa(mb, 'Mass'))
                error('Both inputs must be Mass objects');
            end
            
            mpsa = ma.mps;
            mpsb = mb.mps;
            
            mpsout = [mpsa; mpsb];
            mout = Mass(mpsout);
        end
        
        function [] = plot(mass, varargin)
            N = length(mass.mps);
            x = zeros(N,1);
            y = zeros(N,1);
            z = zeros(N,1);
            for n = 1:N
                pos = mass.mps(n).Pos;
                x(n) = pos(1);
                y(n) = pos(2);
                z(n) = pos(3);
            end
            
            scatter3(x,y,z);
        end
     end
     
     methods (Access = private)
         function [] = computeM(mass)             
             Nps = length(mass.mps);
             
             % mass and cg
             ma = 0;
             mapos = [0 0 0];
             mpMass = zeros(Nps, 1);
             mpPos = zeros(Nps, 3);
             for n = 1:Nps
                 mp = mass.mps(n);
                 mpMass(n) = mp.Mass;
                 mpPos(n,:) = mp.Pos;
                 
                 ma = ma + mpMass(n);
                 mapos = mapos + mpMass(n)*mpPos(n,:);
             end
             
             rnd = 1e8;
             mass.mas = round(rnd*ma)/rnd;
             mass.cg = round(rnd*mapos./ma)./rnd;
             
             if (~isempty(mass.modes))
                 dof = mass.modes.DoF;
                 funcs = mass.modes.MotionFuncs;

                 for n = 1:dof
                     funcs(n).Cg = mass.cg;
                 end
                 
                 pointMotions = cell(dof, Nps);
                 
                 for m = 1:dof
                     for n = 1:Nps
                         pointMotions{m, n} = funcs(m).Evaluate(mpPos(n,:));
                     end
                 end

                 % mass matrix
                 M_ = zeros(dof, dof);
                 for m = 1:dof
                     for n = 1:m
                         iszero = funcs(m).MotionIn*funcs(n).MotionIn';
                         if (iszero ~= 0)
                             mTot = 0;
                             for o = 1:Nps
                                 mval = mpMass(o)*pointMotions{m, o}*pointMotions{n, o}';
                                 mTot = mTot + mval;
                             end
                             M_(m,n) = mTot;
                             M_(n,m) = mTot;
                         end
                     end
                 end

                 rnd = 10^round(log10(rnd/mass.mas));
                 M_ = round(rnd*M_)./rnd;
                 
                 maxM = max(max(M_));
                 minVal = maxM/1e12;
                 
                 for m = 1:dof
                     for n = 1:dof
                         if (abs(M_(m,n)) < minVal)
                             M_(m,n) = 0;
                         end
                     end
                 end
                 
                 mass.M = M_;
             end
             
             mass.iscomp = true;
         end
     end
end