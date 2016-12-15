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
classdef HydroForces
    % Hydrodynamic forces on a floating body
    
    properties (Access = private)
        t;
        nT;
        beta;
        nB;
        a;
        a0;
        ainf;
        b;
        c;
        ex;
        fk;
        dof;
        h;
        rho;
    end

    properties (Dependent)
        T;              % Periods (s)
        Beta;           % Directions (deg)
        A;              % Added mass matrix (dof x dof)
        A0;             % 0 frequency (infinite period) added mass
        Ainf;           % Infinte frequency (0 period) added mass
        B;              % Hydrodynamic Damping matrix (dof x dof)
        C;              % Hydrostatic stiffness matrix (dof x dof)
        Fex;            % Exciting Forces (dof x 1)
        Ffk;            % Froude-Krylov Forces (dof x 1)
        DoF;            % Number of degrees of freedom
        H;              % Water depth
        Rho;            % Fluid density
    end
    
    methods
        function [hf] = HydroForces(t, beta, A, B, C, Ex, varargin)
            hf.nT = length(t);
            hf.nB = length(beta);
            
            hf.t = t;
            hf.beta = beta;
            
            [nt dof1 dof2] = size(A);
            if ((nt ~= hf.nT) || (dof1 ~= dof2))
                error('Invalid added mass matrix');
            end
            
            dofA = dof1;
            
            [nt dof1 dof2] = size(B);
            if ((nt ~= hf.nT) || (dof1 ~= dof2))
                error('Invalid damping matrix');
            end
            
            if (dofA ~= dof1)
                error('HydroForces matrix size mismatch');
            end
            
            [dof1 dof2] = size(C);
            if (dof1 ~= dof2)
                error('Invalid hydrostatic stiffness matrix');
            end
            
            if (dofA ~= dof1)
                error('HydroForces matrix size mismatch');
            end
            
            [nt nb dofe] = size(Ex);
            if ((nt ~= hf.nT) || (nb ~= hf.nB) || (dofe ~= dofA))
                error('Invalid exciting force vector');
            end
            
            if ~isempty(varargin)
                hf.h = varargin{1};
            else
                hf.h = [];
            end
            
            if length(varargin) >= 2
                hf.rho = varargin{2};
            else
                hf.rho = [];
            end
            
            if length(varargin) >= 3
                hf.fk = varargin{3};
            else
                hf.fk = [];
            end
            
            if length(varargin) >= 4
                hf.a0 = varargin{4};
            else
                hf.a0 = [];
            end
            
            if length(varargin) >= 5
                hf.ainf = varargin{5};
            else
                hf.ainf = [];
            end

            hf.a = A;
            hf.b = B;
            hf.c = C;
            hf.ex = Ex;
            hf.dof = dofe;
        end
        
        function [t_] = get.T(hf)
            % Get the wave periods
            t_ = hf.t;
        end
        
        function [bet] = get.Beta(hf)
            % Get the wave headings
            bet = hf.beta;
        end
        
        function [a_] = get.A(hf)
            % Get the added mass
            a_ = hf.a;
        end
        
         function [a_] = get.A0(hf)
            % Get the zero frequency (infintie period) added mass
            a_ = hf.a0;
         end
        
          function [a_] = get.Ainf(hf)
            % Get the infinite frequency (zero period) added mass
            a_ = hf.ainf;
        end
        
        function [b_] = get.B(hf)
            % Get the hydrodynamic damping
            b_ = hf.b;
        end
        
        function [c_] = get.C(hf)
            % Get the hydrostatic stiffness
            c_ = hf.c;
        end
        
        function [e] = get.Fex(hf)
            % Get the hydrodynamic excitation force
            e = hf.ex;
        end
        
        function [e] = get.Ffk(hf)
            % Get the hydrodynamic Froude-Krylov force
            e = hf.fk;
        end
        
        function [dof_] = get.DoF(hf)
            % Get the number of degrees of freedom
            dof_ = hf.dof;
        end
        
        function [h_] = get.H(hf)
            % Get the water depth
            h_ = hf.h;
        end
        
        function [rh] = get.Rho(hf)
            % Get the fluid density
            rh = hf.rho;
        end
    end
end