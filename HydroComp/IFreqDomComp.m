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
classdef IFreqDomComp < IEnergyComp & handle
    % Hydrodynamic frequency domain compuation interface and abstract class.  
    
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
        fext;
        fpto;
        ipto;
        isComp;
        P;
        badInds;
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
        Fext;           % External excitation forces that do not contribute to PTO power absorption
                        % e.g Morison excitation force
        FextPto;        % External excitation forces that contribute to PTO power absorption
        DoF;            % Degrees of freedom
        PTOInds;        % Index matrix of DoFxDoF where a value of 1 indicates the PTO DoF
        Modes;          % String description of all modes of operation
        Pmat;           % Linear constraint matrix;
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
        
        function [val] = get.Fext(hcomp)
            % External excitation forces that do not contribute to PTO power absorption
            % e.g Morison excitation force
            val = hcomp.fext;
        end
        
        function [val] = get.FextPto(hcomp)
            % External excitation forces that contribute to PTO power absorption
            val = hcomp.fpto;
        end
        
        function [dof_] = get.DoF(hcomp)
            % The number of degrees-of-freedom
            dof_ = hcomp.dof;
        end
        
        function [val] = get.PTOInds(hcomp)
            % Index matrix of DoFxDoF where a value of 1 indicates the PTO DoF
                        
            if isempty(hcomp.ipto)
                if ~isempty(hcomp.dpto)
                    hcomp.ipto = hcomp.dpto > 0;
                end
            end
            
            val = hcomp.ipto;
        end
        
        function [P_] = get.Pmat(hcomp)
            % Linear constraint matrix
            P_ = hcomp.P;
        end
                
        function [motions] = Motions(hcomp, varargin)
            % The complex motion amplitudes of the bodies in the array.
            % Optional input includes 'Optimal' which returns the motions
            % required for optimal power absorption of the array
            
            hcomp.computeIfNot();
            
            opts = checkOptions({{'Optimal'}, {'ConstOpt', 1}, {'OrgCoor'}}, varargin);
            
            optm = (opts(1) || opts(2));
            
            orgCoor = opts(3);
            
            omega = 2*pi./hcomp.t;
            
            fex = hcomp.Fex;
            if ~isempty(hcomp.fext)
                fex = fex + hcomp.fext;
            end
            
            if ~isempty(hcomp.fpto)
                fex = fex + hcomp.fpto;
            end
            
            dfreqPto = false;
            if ndims(hcomp.dpto) == 3
                dfreqPto = true;
            end
            
            dfreqPar = false;
            if ndims(hcomp.dpar) == 3
                dfreqPar = true;
            end
            
            kfreq = false;
            if (ndims(hcomp.k) == 3)
                kfreq = true;
            end
            
            if (~optm)
                motions = zeros(size(fex));
                m_ = hcomp.m;
                ddpto = hcomp.dpto;
                ddpar = hcomp.dpar;
                kk = hcomp.k;
                c_ = hcomp.c;
                
                for n = 1:hcomp.nT
                    
                    a_ = reshape(hcomp.a(n,:,:), [hcomp.dof, hcomp.dof]);
                    b_ = reshape(hcomp.b(n,:,:), [hcomp.dof, hcomp.dof]);
                    
                    if (dfreqPto)                       
                        d_ = reshape(ddpto(n,:,:), [hcomp.dof, hcomp.dof]);
                    else
                        d_ = ddpto;
                    end
                    
                    if (dfreqPar)
                        dp = reshape(ddpar(n,:,:), [hcomp.dof, hcomp.dof]);
                        d_ = d_ + dp;
                    else
                        d_ = d_ + ddpar;
                    end
                    
                    if (kfreq)
                        k_ = reshape(kk(n,:,:), [hcomp.dof, hcomp.dof]);
                    else
                        k_ = kk;
                    end

                    lhs = -omega(n)^2.*(m_ + a_) + 1i*omega(n)*(d_ + b_) + k_ + c_;

                    if (ndims(fex) == 2)
                        f = fex(n, :).';
                        if (~hcomp.badInds(n))
                            motions(n, :) = lhs\f;
                        else
                            motions(n, :) = NaN;
                        end
                    else
                        for j = 1:hcomp.nInc
                            f = zeros(hcomp.dof,1);
                            for p = 1:hcomp.dof
                                f(p) = fex(n, j, p);
                            end
                            %f = squeeze(fex(n, j, :));
                            if (~hcomp.badInds(n))
                                motions(n, j, :) = lhs\f;
                            else
                                motions(n, :) = NaN;
                            end
                        end
                    end
                end
                
                if (orgCoor && ~isempty(hcomp.P))
                    mot1 = motions;
                    PT = hcomp.P.';
                    Ndim = size(PT, 1);
                    if (ndims(fex) == 2)
                        motions = zeros(hcomp.nT, Ndim);
                    else
                        motions = zeros(hcomp.nT, hcomp.nInc, Ndim);
                    end
                    
                    for m_ = 1:hcomp.nT
                        for n = 1:hcomp.nInc
                            if (ndims(fex) == 2)
                                motions(m_, :) = PT*mot1(m_, :);
                            else
                                Nm = size(PT,2);
                                mot1s = zeros(size(PT,2), 1);
                                for o = 1:Nm
                                    mot1s(o) = mot1(m_, n, o);
                                end
                                motions(m_, n, :) = PT*mot1s;
                                %motions(m_, n, :) = PT*squeeze(mot1(m_, n, :));
                            end
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
                motions = hcomp.Motions(varargin{:});
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
                        if (~hcomp.badInds(n))
                            v = hcomp.computeOptVel(f, b_, G);
                        else
                            v = NaN;
                        end
                        velocities(n, :) = v;
                    else
                        for j = 1:hcomp.nInc
                            f = squeeze(fex(n, j, :));
                            if (~hcomp.badInds(n))
                                v = hcomp.computeOptVel(f, b_, G);
                            else
                                v = NaN;
                            end
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
            
            [opts, args] = checkOptions({{'Optimal'}, {'ConstOpt', 1}}, varargin);
            
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
                        u = zeros(hcomp.dof, 1);
                        for p = 1:hcomp.dof
                            u(p) = vel(n, j, p);
                        end
                        
                        if ~isempty(hcomp.fpto)
                            fp = zeros(hcomp.dof, 1);
                            for p = 1:hcomp.dof
                                fp(p) = hcomp.fpto(n, j, p);
                            end
                            power(n, j, :) = 0.5*real((d_*conj(u)).*u) + 0.5*real(conj(fp).*u);
                        else
                            power(n, j, :) = 0.5*real((d_*conj(u)).*u);
                        end
                    end
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
        
        function [power] = PowerRAO(hcomp, varargin)      
            % Power response per wave amplitude squared in kW/m^2
            % or if option 'CW' is given, then it is the capture width
            
            [opts, args] = checkOptions({{'CW', 1}, {'spectrum', 1}}, varargin);
            
            compCW = false;
            rho = [];
            if (opts(1))
                compCW = true;
                rho = args{1};
            end
            
            spectrum = [];
            if opts(2)
                spectrum = args{2};
            end

            powFull = hcomp.Power(varargin);
            inds = hcomp.PTOInds;
            inds = (sum(inds) > 0)';
            % TODO: only first direction
            power = sum(powFull(:,1,inds), 3);
            
            if ~isempty(spectrum)
                a_ = spectrum.Amplitudes;
                if isrow(a_) && iscolumn(power)
                    a_ = a_.';
                end
                power = a_.^2.*power;
            end
            
            if ~compCW
                power = power./10^3;
            else
                uEf = IWaves.UnitEnergyFlux(rho, hcomp.t, hcomp.h);
                CW = power./uEf;
                
                power = CW;
            end
        end
        
        function [power] = AveragePower(hcomp, spectrum, varargin)
            % Average power of the body in the spectrum in kW
            % or if option 'CW' is given, then it is given in capture width
            [opts, args] = checkOptions({{'CW', 1}}, varargin);
            
            compCW = false;
            rho = [];
            if (opts(1))
                compCW = true;
                rho = args{1};
            end
            
            powRao = hcomp.PowerRAO;
            
            amps = spectrum.Amplitudes;
            if (iscolumn(powRao) && isrow(amps)) ...
                    || (isrow(powRao) && iscolumn(amps))
                amps = amps.';
            end
            
            pow0 = powRao.*amps.^2;
            
            power = sum(sum(pow0));
            
            if compCW
                Ef = spectrum.EnergyFlux(rho, hcomp.h)./1000;  % energy flux in kW
                power = power/sum(sum(Ef));
            end
        end
        
        function [energy, energyMat, powerMat, idptos] = AnnualEnergyProd(hcomp, waveClim, varargin)
            
            [opts, args] = checkOptions({{'Sample', 1}, {'ratedPow', 1}, {'dptos', 1}, {'minOcc', 1}}, varargin);
            
            sampSize = 1;
            if opts(1)
                sampSize = args{1};
            end
            
            perOccur = false;
            ratedPow = [];
            if opts(2)
                ratedPow = args{2};
                perOccur = true;
            end
            
            dptos = [];
            if opts(3)
                dptos = args{3};
                perOccur = true;
            end
            hrsYr = 24*365;
            occlim = 1/hrsYr;     % ignore sea states that occur less than an hour per year
            if opts(4)
                occlim = args{4};
            end
            
            idptos = 1;
            
            if ~perOccur
                avgSpec = waveClim.AverageSpectrum;
                % average power in kW
                pow = hcomp.AveragePower(avgSpec);

                % MWh/yr
                energy = pow*hrsYr./10^3;
            
            else
                plim = 1e3;             % ignore sea states below a 1 kW total
                
                [powerMat, idptos] = PowerMatrix.CreatePowerMatrix(hcomp, [], [], ...
                    'waveClim', waveClim, 'minPow', plim, 'minOcc', occlim, ...
                    'makeObj', 'ratedPow', ratedPow, 'dptos', dptos);
                
                Pow = powerMat.Matrix;
                
                freqOccs = waveClim.FreqOccurance;
                energyMat = hrsYr*freqOccs.*Pow; % kWh/yr

                % Total power: MWh/yr
                energy = sum(sum(energyMat))./1e3;
            end
            
            energy = energy*ones(sampSize);
        end
        
        function [pmat, idptos] = PowerMatrix(hcomp, Hs, T, varargin)
            
            [pmat, idptos] = PowerMatrix.CreatePowerMatrix(hcomp, Hs, T, 'makeObj', varargin{:});
        end
        
        function [frad] = Frad(hcomp, varargin)
            
            opts = checkOptions({{'Optimal'}, {'ConstOpt', 1}, {'OrgCoor'}}, varargin);
            
            if opts(3)
                xi = hcomp.Motions;
                orgCoor = true;
            else
                xi = hcomp.Motions(varargin{:});
                orgCoor = false;
            end
            omega = 2*pi./hcomp.t;
            
            frad = zeros(size(xi));
            
            for n = 1:hcomp.nT
                a_ = zeros(hcomp.dof);
                b_ = zeros(hcomp.dof);
                
                for p = 1:hcomp.dof
                    for q = 1:hcomp.dof
                        a_(p,q) = hcomp.a(n,p,q);
                        b_(p,q) = hcomp.b(n,p,q);
                    end
                end
                if ndims(xi) == 2
                    frad(n,:) = (-(-omega(n)^2.*a_ + 1i*omega(n)*b_)*xi(n,:).').';
                else
                    for j = 1:hcomp.nInc
                        frad(n,j,:) = (-(-omega(n)^2.*a_ + 1i*omega(n)*b_)*squeeze(xi(n,j,:)));
                    end
                end
            end
            
            if (orgCoor && ~isempty(hcomp.P))
                frad1 = frad;
                PT = hcomp.P.';
                Ndim = size(PT, 1);
                if (ndims(frad1) == 2)
                    frad = zeros(hcomp.nT, Ndim);
                else
                    frad = zeros(hcomp.nT, hcomp.nInc, Ndim);
                end

                for m_ = 1:hcomp.nT
                    for n = 1:hcomp.nInc
                        if (ndims(frad1) == 2)
                            frad(m_, :) = PT*frad1(m_, :);
                        else
                            frad(m_, n, :) = PT*squeeze(frad1(m_, n, :));
                        end
                    end
                end
            end
        end
        
        function [] = SetM(hcomp, m_)
            % Set the Mass value, currently does not change the
            % values of the floating bodies
            
            hcomp.checkMatSize(m_);

            hcomp.m = m_;
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
        
        function [] = SetPTOInds(hcomp, val)
            % Set the index matrix of DoFxDoF where a value of 1 indicates the PTO DoF
            
%             if ndims(val) > 2
%                 error('PTOInds must be a Ndof x 1 vector');
%             end
%             
%             [Nx, Ny] = size(val);
%             if Nx ~= hcomp.dof || Ny ~= 1
%                  error('PTOInds must be a Ndof x 1 vector');
%             end
%             
%             if (sum(val == 0) + sum(val == 1)) ~= (hcomp.dof)
%                 error('The PTOInds matrix must be 1''s or 0''s');
%             end
            
            hcomp.ipto = val;
        end
        
        function [] = SetFext(hcomp, val)
            % External excitation forces that do not contribute to PTO power absorption
            % e.g Morison excitation force
            if ~isempty(val)
                if ~all(size(hcomp.Fex) == size(val))
                    error('Fext must be the same size as Fex');
                end
            end
            hcomp.fext = val;
        end
        
        function [] = SetFpto(hcomp, val)
            % External excitation forces that contribute to PTO power absorption
            if ~isempty(val)
                if ~all(size(hcomp.Fex) == size(val))
                    error('Fpto must be the same size as Fex');
                end
            end
            hcomp.fpto = val;
        end
    end
    
    methods (Access = protected)        
        function [] = initHydroParam(hcomp, t_, h_, bods, const, P)    
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
            
            dof_ = IFreqDomComp.GetDoF(bods);
            if (const)
                [M, N] =  size(P);
                if(N ~= dof_)
                    error('Constraint P matrix not correct size');
                end
                dof_ = M;
                hcomp.P = P;
            end
            
            
            [m_, dpto_, dpar_, k_, c_] = IFreqDomComp.resizeMDK(bods);

            if (const)
                m_ = P*m_*P.';
                dpto_ = P*dpto_*P.';
                dpar_ = P*dpar_*P.';
                k_ = P*k_*P.';
                c_ = P*c_*P.';
            end
            hcomp.m = m_;
            hcomp.dpto = dpto_;
            hcomp.dpar = dpar_;
            hcomp.k = k_;
            hcomp.c = c_;

            hcomp.dof = dof_;
            
            % External excitation force;
            hcomp.fext = [];
            hcomp.fpto = [];
            
            hcomp.badInds = zeros(hcomp.nT);
            
            % PTO indices
            hcomp.ipto = [];
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
             %v = 0.5*inv(b_)*f;
              v = 0.5*(b_\f);
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
        
        function [] = checkBadVals(hcomp, B)            
            for m_ = 1:hcomp.nT
                b_ = squeeze(B(m_,:,:));
                
                [N, ~] = size(b_);
                
                for n = 1:N
                    if (b_(n,n) < 0)
                        hcomp.badInds(m_) = 1;
                    end
                end
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
            
            if (isempty(geo.C))
                evalC = false;
                 c_ = zeros(df, df);
            else
                cdof = size(geo.C,1);
                if (cdof == df)
                    evalC = false;
                    c_ = geo.C;
                else
                    evalC = true;
                    c_ = zeros(df, df);
                end
            end

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
%                         if (evalC)
%                             c_(lsf + j, lsf + p) = geo.C(iv(j), iv(p));
%                         end
                    end
                end
                if evalC
                    cdof = size(geo.C,1);
                    for j = 1:cdof
                        for p = 1:cdof
                            c_(lsf + j, lsf + p) = geo.C(j, p);
                        end
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