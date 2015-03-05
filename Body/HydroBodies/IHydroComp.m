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
classdef IHydroComp < handle
    % Hydrodynamic compuation interface and abstract class.  
    
    properties (SetAccess = protected, GetAccess = protected)
        t;
        nT;
        iwaves;
        nInc;
        h;
        dof;
        a;
        b;
        c;
        m;
        dpto;
        dpar;
        k;
        kgm;
        isComp;
    end
    
    properties (Dependent)
        T;              % Periods (s)
        H;              % Water depth
        A;              % Added mass matrix (dof x dof)
        B;              % Hydrodynamic Damping matrix (dof x dof)        
        C;              % Hydrostatic stiffness matrix (dof x dof)
        M;              % Mass matrix for all bodies
        Dpto;           % PTO Damping matrix for all bodies 
        Dpar;           % Parasitic damping matrix for all bodies
        K;              % Mechanical Stiffness matrix for all bodies
        DoF;            % Degrees of freedom
        Modes;          % String description of all modes of operation
    end
    
    properties (Abstract)
        IncWaves;       % Incident Waves
        Bodies;         % Individual bodies in array
        Fex;            % Exciting Forces 
    end
    
    methods (Abstract, Access = protected)   
        computeIfNot(hcomp);
    end
    
    methods        
        function [t_] = get.T(hcomp)
            % The wave periods
            t_ = hcomp.t;
        end
                
        function [h_] = get.H(hcomp)
            % The water depth
            h_ = hcomp.h;
        end
        
        function [a_] = get.A(hcomp)
            % The hydrodynamic added mass
            hcomp.computeIfNot();
            
            a_ = hcomp.a;
        end
        
        function [b_] = get.B(hcomp)
            % The hydrodynamic damping
            hcomp.computeIfNot();
            
            b_ = hcomp.b;
        end
        
        function [c_] = get.C(hcomp)
            % The hydrostatic stiffness
            hcomp.computeIfNot();
            
            c_ = hcomp.c;
        end
                
        function [m_] = get.M(hcomp)
            % The mass matrix
            m_ = hcomp.m;
        end
        
        function [d_] = get.Dpto(hcomp)
            % The power take-off damping
            d_ = hcomp.dpto;
        end
        
        function [d_] = get.Dpar(hcomp)
            % The parastic damping
            d_ = hcomp.dpar;
        end
        
        function [k_] = get.K(hcomp)
            % The mechanical stiffness
            k_ = hcomp.k;
        end
        
        function [dof_] = get.DoF(hcomp)
            % The number of degrees-of-freedom
            dof_ = hcomp.dof;
        end
                
        function [motions] = Motions(hcomp, varargin)
            % The complex motion amplitudes of the bodies in the array.
            % Optional input includes 'Optimal' which returns the motions
            % required for optimal power absorption of the array
            
            hcomp.computeIfNot();
            
            opts = checkOptions({{'Optimal'}, {'ConstOpt', 1}}, varargin);
            
            optm = (opts(1) || opts(2));
            
            omega = 2*pi./hcomp.t;
            
            fex = hcomp.Fex;
            
            dfreq = false;
            if (ndims(hcomp.dpto) == 3)
                dfreq = true;
            end
            
            kfreq = false;
            if (ndims(hcomp.k) == 3)
                kfreq = true;
            end
            
            if (~optm)
                motions = zeros(size(fex));
                m_ = hcomp.m;
                dd = hcomp.dpto + hcomp.dpar;
                kk = hcomp.k;
                c_ = hcomp.c;
                
                for n = 1:hcomp.nT
                    a_ = squeeze(hcomp.a(n,:,:));
                    b_ = squeeze(hcomp.b(n,:,:));
                    
                    if (dfreq)
                        d_ = squeeze(dd(n,:,:));
                    else
                        d_ = dd;
                    end
                    
                    if (kfreq)
                        k_ = squeeze(kk(n,:,:));
                    else
                        k_ = kk;
                    end

                    lhs = -omega(n)^2.*(m_ + a_) + 1i*omega(n)*(d_ + b_) + k_ + c_;

                    if (ndims(fex) == 2)
                        f = fex(n, :).';
                        motions(n, :) = lhs\f;
                    else
                        for j = 1:hcomp.nInc
                            f = squeeze(fex(n, j, :));
                            motions(n, j, :) = lhs\f;
                        end
                    end
                end
            else
                vel = hcomp.Velocities(varargin{:});
                motions = zeros(size(vel));
                
                for n = 1:hcomp.nT
                    if (ndims(fex) == 2)
                        motions(n,:) = -1i/omega(n)*vel(n,:);
                    else
                        motions(n,:,:) = -1i/omega(n)*vel(n,:,:);
                    end
                end
            end
        end
        
        function [velocities] = Velocities(hcomp, varargin)
            % The complex velocity amplitudes of the bodies in the array
            % Optional input includes 'Optimal' which returns the
            % velocities required for optimal power absorption of the 
            % array
            
            hcomp.computeIfNot();
            
            [opts, args] = checkOptions({{'Optimal'}, {'ConstOpt', 1}}, varargin);
            
            if (opts(1))
                optm = true;
            else
                optm = false;
            end
            
            if (opts(2))
                coptm = true;
                gamma = args{2};
            else
                coptm = false;
            end
            
            if (optm && coptm)
                error('Either ''Optimal'' or ''ConstOpt'' may be used as an argument, not both');
            end
            
            if (coptm)
                nGam = length(gamma);
                if (nGam ~= hcomp.dof)
                    error('The constaint vector must be the same size as the DOF');
                end
                
                G = diag(gamma);
                optm = true;
            else
                G = NaN;
            end
            
            fex = hcomp.Fex;
            
            if (~optm)
                motions = hcomp.Motions;
                omega = 2*pi./hcomp.t;

                velocities = zeros(size(motions));

                for n = 1:hcomp.nT
                    if (ndims(fex) == 2)
                        velocities(n,:) = 1i*omega(n)*motions(n,:);
                    else
                        velocities(n,:,:) = 1i*omega(n)*motions(n,:,:);
                    end
                end
            else
                % optimal motions
                % or
                % constained optimal velocity following 
                % Pizer (1993) "Maximum wave-power absorption of point
                % absorbers under motion constaints"
                velocities = zeros(size(fex));
                
                for n = 1:hcomp.nT
                    b_ = squeeze(hcomp.b(n,:,:));

                    if (ndims(fex) == 2)
                        f = fex(n, :).';
                        v = hcomp.computeOptVel(f, b_, G);
                        velocities(n, :) = v;
                    else
                        for j = 1:hcomp.nInc
                            f = squeeze(fex(n, j, :));
                            v = hcomp.computeOptVel(f, b_, G);
                            velocities(n, j, :) = v;
                        end
                    end
                end                
            end
        end
        
        function [power] = Power(hcomp, varargin)
            % The power produces by each mode of motion
            % Optional input includes 'Optimal' which returns the
            % velocities required for optimal power absorption of the 
            % array
            
            hcomp.computeIfNot();
            
            opts = checkOptions({{'Optimal'}, {'ConstOpt', 1}}, varargin);
            
            optm = (opts(1) || opts(2));
                        
            if (~optm)
                vel = hcomp.Velocities;

                power = zeros(size(vel));
                
                dfreq = false;
                if (ndims(hcomp.dpto) == 3)
                    dfreq = true;
                end

                for n = 1:hcomp.nT    
                    if (dfreq)
                        % Only the PTO damping is used to compute power
                        d_ = squeeze(hcomp.dpto(n,:,:));
                    else
                        d_ = hcomp.dpto;
                    end
                    for j = 1:hcomp.nInc
                        u = squeeze(vel(n, j, :));
                        power(n, j, :) = 0.5*real((d_*conj(u)).*u);
                    end
%                     if (hcomp.nB == 1)
%                         u = squeeze(vel(n,1,:));
%                         power(n, :) = 0.5*real((hcomp.d*conj(u)).*u);
%                     else
%                         for j = 1:hcomp.nB
%                             u = squeeze(vel(n, j, :));
%                             power(n, j, :) = 0.5*real((hcomp.d*conj(u)).*u);
%                         end
%                     end
                end
            else
                vel = hcomp.Velocities(varargin{:});
                vel = squeeze(vel);
                
                power = zeros(size(vel));
                
                for n = 1:hcomp.nT
                    b_ = squeeze(hcomp.b(n,:,:));

                    if (hcomp.nInc == 1)
                        u = squeeze(vel(n,:)).';
                        power(n, :) = 0.5*real((b_*conj(u)).*u);
                    else
                        for j = 1:hcomp.nInc
                            u = squeeze(vel(n, j, :));
                            power(n, j, :) = 0.5*real((b_*conj(u)).*u);
                        end
                    end
                end
            end
        end
        
        function [] = SetDpto(hcomp, d_)
            % Set the PTO damping values, currently does not change the
            % values of the floating bodies
            
            hcomp.checkMatSize(d_);

            hcomp.dpto = d_;
        end
        
        function [] = SetDpar(hcomp, d_)
            % Set the parasitic damping values, currently does not change the
            % values of the floating bodies
            
            hcomp.checkMatSize(d_);
            
            hcomp.dpar = d_;
        end
        
        function [] = SetK(hcomp, k_)
            % Set the mechanical stiffness values, currently does not change the
            % values of the floating bodies
            
            hcomp.checkMatSize(k_);
            
            hcomp.k = k_;
        end
    end
    
    methods (Access = protected)        
        function [] = initHydroParam(hcomp, t_, h_, bods)    
            if (isvector(t_))
                if (isrow(t_))
                    t_ = t_.';
                end
            else
                error('Periods must be a vector');
            end
            hcomp.t = t_;
            hcomp.nT = length(hcomp.t);
            
            hcomp.h = h_;
            
            [m_, dpto_, dpar_, k_, c_] = IHydroComp.resizeMDK(bods);

            hcomp.m = m_;
            hcomp.dpto = dpto_;
            hcomp.dpar = dpar_;
            hcomp.k = k_;
            hcomp.c = c_;

            hcomp.dof = IHydroComp.GetDoF(bods);
        end
        
        function [] = setIncWaves(hcomp, iwavs)
            nIwav = length(iwavs);
            
            for n = 1:nIwav
                if (~isa(iwavs(n), 'IWaves'))
                    error('Incident wave must be an IWaves');
                end

                if(~iwavs(n).IsIncident)
                    error('Waves must be incident waves');
                end

                if any(abs(iwavs(n).T - hcomp.t) > 1e-9)
                    error('Incident waves must have the same wave periods as the HydroComp');
                end

                if (abs(iwavs(n).H - hcomp.h) > 1e-9)
                    error('Incident waves must have the same water depth as the HydroComp');
                end
            end
            
            hcomp.iwaves = iwavs;
            hcomp.nInc = nIwav;
            hcomp.isComp = false;
        end
        
        function [v] = computeOptVel(hcomp, f, b_, G)
             % optimal motions
             % or
             % constained optimal velocity following 
             % Pizer (1993) "Maximum wave-power absorption of point
             % absorbers under motion constaints"
             
             % The constrained optimization is still a bit choppy - I think
             % it is due to numerical errors in computing the eigenvalues.
             
             % This is the unconstained optimal velocity
             v = 0.5*inv(b_)*f;
             if (~isnan(G(1,1)))
                 % want to compute the constained optimal velocity
                 % first, check to see whether the velocity
                 % violates the constaint
                 res = v'*inv(G).^2*v;
                 if (res > 1)
                     % it violate constaint, must compute new velocity
                     GBG = G*b_*G;
                     [Q, Lam] = eig(GBG);
                     Xpri = Q'*G*f;
                     I = eye(hcomp.dof);
                     
                     fun1 = @(x) hcomp.fmu(Xpri, Lam, x);
                     fun = @(x) real(Xpri'*inv(Lam + x*I).^2*Xpri - 4);
                     
                     fstart = 0;
                     fend = 1e6;
                     
                     fun(fstart)
                     fun(fend)
                     
                     fun1(fstart)
                     fun1(fend)
                     
                     mu = fzero(fun, [0 1e6]);
                     
                     v = 0.5*inv(b_ + mu*I)*f;
                 end 
             end    
        end
        
        function [f] = fmu(hcomp, Xpri, Lam, mu)
            N = length(Xpri);
            f = 0;
            for n = 1:N
                f = f + abs(Xpri(n)).^2/(Lam(n,n) + mu).^2;
            end
            f = f- 4;
            f = real(f); % this should always be real, but sometimes it has a very small imag
        end
        
        function [] = checkMatSize(hcomp, m)
            if (ndims(m) == 2)
                [row, col] = size(m);
                if (row ~= hcomp.dof || col ~= hcomp.dof)
                    error('The matrix must be of size DoF x DoF');
                end
            elseif (ndims(m) == 3)
                [nt, row, col] = size(m);
                if (nt ~= hcomp.nT)
                    error('The number of matrices must be equal to the number of periods');
                end
                if (row ~= hcomp.dof || col ~= hcomp.dof)
                    error('The matrix must be of size DoF x DoF');
                end
            else
                error('Matrix wrong size');
            end
        end
    end
       
    methods (Static)
        function [dof] = GetDoF(fbs)
            % Computes the degrees of freedom from a vector of floating
            % bodies
            nbody = length(fbs);
            dof = 0;

            for n = 1:nbody
                fb = fbs(n);
                dof = dof + fb.Modes.DoF;
            end
        end
        
        function [a] = IAmps(M, pos, k0, beta)
            
            A0 = exp(-1i*k0*(pos(1)*cos(beta) + pos(2)*sin(beta)));

            a = zeros(2*M+1,1);

            for m = -M:M
                a(m+M+1) = A0*exp(-1i*m*(beta + pi/2));
            end
        end
    end
    
    methods (Static, Access = protected)
                
        function [m_, dpto_, dpar_, k_, c_] = resizeMDK(fbs)
            nbody = length(fbs);
            % Get Mass, Damping and Stiffness matrices from geomerties.
            % Can't handle connected bodies...
            df = 0;

            for n = 1:nbody
                geo = fbs(n);
                df = df + geo.Modes.DoF;
            end

            m_ = zeros(df, df);
            dpto_ = zeros(df, df);
            dpar_ = zeros(df, df);
            k_ = zeros(df, df);
            c_ = zeros(df, df);

            lsf = 0;
            for n = 1:nbody
                geo = fbs(n);
                v = geo.Modes.Vector;
                count = sum(v);
                iv = find(v == 1);
                for j = 1:count
                    for p = 1:count
                        m_(lsf + j, lsf + p) = geo.M(iv(j), iv(p));
                        dpto_(lsf + j, lsf + p) = geo.Dpto(iv(j), iv(p));
                        dpar_(lsf + j, lsf + p) = geo.Dpar(iv(j), iv(p));
                        k_(lsf + j, lsf + p) = geo.K(iv(j), iv(p));
                        c_(lsf + j, lsf + p) = geo.C(iv(j), iv(p));
                    end
                end
                lsf = lsf + count;
            end
        end

        function [modes] = getModes(fbs)        
            % Computes string of active modes of floating bodies
            df = FloatingBodyArray.GetDoF(fbs);
            modes = cell(df, 1);
            
            nbody = length(fbs);

            fbNames = cell(nbody, 1);

            for n = 1: nbody
                fb = fbs(n);
                fbName = fb.Handle;
                for o = 1:(n-1)
                    if (strcmp(fbNames{o}, fbName))
                        error('All floating bodies in a given array must have unique handles.');
                    end
                end
                fbNames{n} = fbName;
                
                mo = fb.Modes.Motions;
                cnt = length(mo);
                
                for o = 1:cnt;
                    modes{n + o - 1} = [fbName ' - ' mo{o}];
                end
            end
        end
    end
end