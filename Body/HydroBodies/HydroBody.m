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
classdef HydroBody < FloatingBody
    % A floating body with attached hydrodynamic forces - added mass, 
    % damping, excitation force, a diffraction transfer matrix, a force 
    % excitation matrix etc. 
    %
    % Shall be used to compute hydrodnamic array interactions forces.
    
    properties (Access = protected)
        nT;
        hydroForces;
        diffTM;
        fexTM;
        radCoefs;
        mlim;
    end

    properties (Dependent)
        T;              % Periods (s)
        A;              % Added mass matrix (dof x dof)
        B;              % Hydrodynamic Damping matrix (dof x dof)
        DoF;            % Number of degrees of freedom
        H;              % The water depth
        Rho;            % The fluid density
        DiffTM;         % Diffraction transfer matrix
        FexTM;          % Force transfer matrix
        RadCoefs;       % Radiation coefficients
        Mlim;           % The M limit value for the partial circular wave summation
    end
 
    methods

        function [hb] = HydroBody(body, varargin)
            % Constructor
            if (~isa(body,'FloatingBody'))
                error('fb must be a FloatingBody object');
            end
            
            hb = hb@FloatingBody(body);
            
            if (isa(body, 'HydroBody'))
                hb.nT = length(body.hydroForces.T);
                hb.hydroForces = body.hydroForces;
                hb.diffTM = body.diffTM;
                hb.fexTM = body.fexTM;
                hb.radCoefs = body.radCoefs;
                hb.mlim = body.mlim;
            else            
                if (length(varargin) ~= 4)
                    error('If not using the copy constructor, there must be 5 arguments: body, hydF, dtm, ftm, radCs')
                end
                
                hydF = varargin{1};
                dtm = varargin{2};
                ftm = varargin{3};
                radCs = varargin{4};

                if (~isa(hydF, 'HydroForces'))
                    error('hydF must be a HydroForces object');
                end

                hb.hydroForces = hydF;
                hb.nT = length(hydF.T);

                % One DTM and FTM for each wave period
                dof = hydF.DoF;
                if (length(dtm) ~= hb.nT)
                    error('The number of diffraction transfer matrices must be the same as the number of wave periods.');
                end

                if (length(ftm) ~= hb.nT)
                    error('The number of force transfer matrices must be the same as the number of wave periods.');
                end
                
                % One set of radiation coefficients fo each DoF
                [Nt, Q] = size(radCs);
                if ((Nt ~= hb.nT) || (Q ~= hb.hydroForces.DoF))
                    error('The number of sets of radiation coefficients must be equal to the number periods by the number of degrees of freedom.');
                end

                hb.mlim = zeros(hb.nT, 1);
                
                for m = 1:hb.nT
                    % Each DTM must be square
                    [M, N] = size(dtm{m});
                    if (M ~= N)
                        error('All diffraction transfer matrices must be square');
                    end

                    % The FTM must be DoF x 2M+1
                    [Q, M] = size(ftm{m});
                    if (Q ~= dof)
                        error('All force transfer matrices must be DoF x 2N+1');
                    end
                    
                    if (M ~= N)
                        error('The diffraction transfer matix and the force transfer matrix must have the same M limits');
                    end
                    
                    % Rad Coefs
                    for n = 1:hb.DoF
                        if (size(radCs{m,n}) ~= [1, M])
                            error('The radiation coefficient must have the same M limits as the diffraction transfer and force transfer matrices');
                        end
                    end
                    
                    hb.mlim(m) = (M-1)/2;
                end

                hb.diffTM = dtm;
                hb.fexTM = ftm;
                hb.radCoefs = radCs;
            end
        end
        
        function [t_] = get.T(hb)
            % Get the wave periods
            t_ = hb.hydroForces.T;
        end
        
        function [a_] = get.A(hb)
            % Get the added mass
            a_ = hb.hydroForces.A;
        end
        
        function [b_] = get.B(hb)
            % Get the hydrodynamic damping
            b_ = hb.hydroForces.B;
        end
                
        function [dof] = get.DoF(hb)
            % Get the number of degrees of freedom
            dof = hb.hydroForces.DoF;
        end
        
        function [h_] = get.H(hb)
            % Get the water depth
            h_ = hb.hydroForces.H;
        end
        
        function [rh] = get.Rho(hb)
            % Get the fluid density
            rh = hb.hydroForces.Rho;
        end
        
        function [dtm] = get.DiffTM(hb)
            % Get the diffraction transfer matrix
            if (hb.angle == 0)
                dtm = hb.diffTM;
            else
                dtm = cell(hb.nT, 1);
                psi = pi/180*hb.angle;
                for n = 1:hb.nT
                    dtmn = hb.diffTM{n};
                    [N ~] = size(dtmn);
                    dtmn2 = zeros(N,N);
                    N = (N-1)/2;
                    for p = -N:N
                        for q = -N:N
                            dtmn2(p+N+1, q+N+1) = exp(1i*(q-p)*psi)*dtmn(p+N+1, q+N+1);
                        end
                    end
                    dtm{n} = dtmn2;
                end
            end
        end
        
        function [ftm] = get.FexTM(hb)
            % Get the force transfer matrix
            if (hb.angle == 0)
                ftm = hb.fexTM;
            else
                ftm = cell(hb.nT, 1);
                psi = pi/180*hb.angle;
                for n = 1:hb.nT
                    ftmn = hb.fexTM{n};
                    [N M] = size(ftmn);
                    ftmn2 = zeros(N,M);
                    N = (M-1)/2;
                    for p = 1:hb.DoF
                        for q = -N:N
                            ftmn2(p, q+N+1) = exp(1i*q*psi)*ftmn(p, q+N+1);
                        end
                    end
                    ftm{n} = ftmn2;
                end
            end
        end
        
        function [rcfs] = get.RadCoefs(hb)
            % Get the radiation coefficients
            if (hb.angle == 0)
                rcfs = hb.radCoefs;
            else
                dof = hb.hydroForces.DoF;
                rcfs = cell(hb.nT, dof);
                psi = pi/180*hb.angle;
                for n = 1:hb.nT
                    for m = 1:dof
                        rcfnm = hb.radCoefs{n, m};
                        [M, N] = size(rcfnm);
                        rcfnm2 = zeros(1, N);
                        N = (N-1)/2;
                        for p = -N:N
                            rcfnm2(p+N+1) = exp(-1i*p*psi)*rcfnm(p+N+1);
                        end
                     
                        rcfs{n, m} = rcfnm2;
                    end
                end
            end
        end
        
        function [ml] = get.Mlim(hb)
            ml = hb.mlim;
        end
        
        function [areEq] = eq(hba, hbb)
            % Overloaded equality operator
            areEq = 1;
            
            if (~isa(hba, 'HydroBody') || ~isa(hbb, 'HydroBody'))
                error('Each argument must be a Hydrobody');
            end
            
            % Check periods
            if any(hba.T ~= hbb.T)
                areEq = 0;
                return;
            end
                        
            % Check depth
            if ~(isempty(hba.H) && isempty(hbb.H))
                if (hba.H ~= hbb.H)
                    areEq = 0;
                    return;
                end
            end
            
            % Check density
            if ~(isempty(hba.Rho) && isempty(hbb.Rho))
                if (hba.Rho ~= hbb.Rho)
                    areEq = 0;
                    return;
                end
            end
        end
        
        function [areNE] = ne(hba, hbb)
            % Overloaded inequality operator
            areNE = 1;
            if (hba == hbb)
                areNE = 0;
            end
        end
    end
    
    methods (Static)
        
        function [MatOut] = Resize2M(MatIn, M)
            [nrow, ncol] = size(MatIn);
            
            Min = ncol;
            if (mod(Min,2) ~= 1)
                error('The number of cols of the input matrix must be odd');
            end
            Min = (Min - 1)/2;
            
            if (nrow == 1)
                if (M > Min)
                    MatOut = zeros(1, 2*M+1);
                    MatOut((M + 1 - Min):(M + 1 + Min)) = MatIn;
                elseif (M < Min)
                    MatOut = MatIn((Min + 1 - M):(Min + 1 + M));
                else
                    MatOut = MatIn;
                end
            elseif (nrow == ncol)
                if (M > Min)
                    MatOut = zeros(2*M+1,2*M+1);
                    MatOut((M + 1 - Min):(M + 1 + Min),(M + 1 - Min):(M + 1 + Min)) = MatIn;
                elseif (M < Min)
                    MatOut = MatIn((Min + 1 - M):(Min + 1 + M), (Min + 1 - M):(Min + 1 + M));
                else
                    MatOut = MatIn;
                end
            else
                error('Input matrix not of recognizable size');
            end
        end
    end
end