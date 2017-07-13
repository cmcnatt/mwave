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
classdef FreqDomComp < IFreqDomComp & CopyLoadHandleBase
    % Computes properties (motions, velocity, power) for a single or an 
    % array of floating bodies in the frequency domain
    
    properties (Access = private)
        fbs;
        isHB;
        fex;
        ffk;
        devCount;
    end
    
    properties (Dependent)
        IncWaves;
        Bodies;
        Fex;
        Ffk;
        DeviceCount;
    end
    
    methods
        
        function [hbcomp] = FreqDomComp(varargin)
            % Constructor - takes either a single HydroBody of a FreqDomForces object and and a vector of
            % floating bodies
            % If arg1 is a HydroBody, then arg 2 is IWaves
            % If arg1 is a FreqDomForces, then arg 2 is Floating Body
            
            [opts, args] = checkOptions({{'UseBodC'}, {'Constrained', 1}, {'CheckNegB'}}, varargin);
            useBodC = opts(1);
            const = opts(2);
            if (const)
                P = args{2};
            else
                P = [];
            end
            checkNegB = opts(3);
            
            hbcomp.isHB = false;
            
            if (~isempty(varargin))
                hydro = varargin{1};
                
                if (isa(hydro, 'HydroBody'))
                    hbcomp.isHB = true;
                    if (isempty(varargin{2}))
                        error('Second argument must be IWaves');
                    elseif(~isa(varargin{2}, 'IWaves'));
                        error('Second argument must be IWaves');
                    end
                    if (const)
                        error('Cannot created a constrained computation with a HydroBody');
                    end
                    iwavs = varargin{2};
                    useBodC = true;
                elseif (isa(hydro, 'FreqDomForces'))                
                    if (isempty(varargin{2}))
                        error('Second input must be a FloatingBody');
                    elseif(~isa(varargin{2}, 'FloatingBody'));
                        error('Second argument must be FloatingBody');
                    end
                    fbods = varargin{2};
                    
                    if (useBodC && isempty(fbods(1).C))
                        error('Floating body does not have a hydrostatic matrix');
                    end

                    nB = length(hydro.Beta);
                    wuns = ones(size(hydro.T));
                    for n = 1:nB
                        iwavs(n) = PlaneWaves(wuns, hydro.T, hydro.Beta(n)*wuns, hydro.H);
                    end
                else
                    error('Inputs must be either a HydroBody or a FreqDomForces object and and a vector of floating bodies');
                end
                                           
                if (~hbcomp.isHB)       
                    for n = 1:length(fbods)
                        if (~isa(fbods(n), 'FloatingBody'))
                            error('Second input must be a FloatingBody')
                        end
                    end

                    hbcomp.initHydroParam(hydro.T, hydro.H, fbods, const, P);
                    hbcomp.setIncWaves(iwavs);
                    fex_ = hydro.Fex;
                    ffk_ = hydro.Ffk;
                    if (const)
                        [M, N] =  size(P);
                        [Nt, Nb, ~] = size(fex_);
                        fexc = zeros(Nt, Nb, M);
                        if ~isempty(ffk_)
                            ffkc = zeros(Nt, Nb, M);
                        end
                        for m = 1:Nt
                            for n = 1:Nb
                                f = zeros(N, 1);
                                for p = 1:N
                                    f(p) = fex_(m,n,p);
                                end
                                %fexc(m,n,:) = P*squeeze(fex_(m,n,:));
                                fexc(m,n,:) = P*f;
                                if ~isempty(ffk_)
                                    f = zeros(N, 1);
                                    for p = 1:N
                                        f(p) = ffk_(m,n,p);
                                    end
                                    %fexc(m,n,:) = P*squeeze(fex_(m,n,:));
                                    ffkc(m,n,:) = P*f;
                                end
                            end
                        end
                        fex_ = fexc;   
                        if ~isempty(ffk_)
                            ffk_ = ffkc;
                        end
                    end
                    hbcomp.fex = fex_;
                    hbcomp.ffk = ffk_;
                else
                    hbcomp.initHydroParam(hydro.T, hydro.H, hydro, const, P);
                    hbcomp.setIncWaves(iwavs);
                end

                if (~useBodC)
                    c = hydro.C;
                    if (const)
                        c = P*c*P.';
                    end
                    hbcomp.c = c;
                end
                
                a_ = hydro.A;
                b_ = hydro.B;
                
                if checkNegB
                    hbcomp.checkBadVals(b_);
                end
                
                if (const)
                    ac = zeros(Nt, M, M);
                    bc = zeros(Nt, M, M);
                    for m = 1:Nt
                        am = zeros(N);
                        bm = zeros(N);
                        for p = 1:N
                            for q = 1:N
                                am(p,q) = a_(m,p,q);
                                bm(p,q) = b_(m,p,q);
                            end
                        end
%                         ac(m,:,:) = P*squeeze(a_(m,:,:))*P.';
%                         bc(m,:,:) = P*squeeze(b_(m,:,:))*P.';
                        ac(m,:,:) = P*am*P.';
                        bc(m,:,:) = P*bm*P.';
                    end
                    a_ = ac;
                    b_ = bc;
                end
                
                hbcomp.a = a_;
                hbcomp.b = b_;

                %[row, col] = size(squeeze(hbcomp.a(1,:,:)));

                if (hbcomp.dof ~= size(hbcomp.m, 1) && hbcomp.dof ~= size(hbcomp.m, 2))
                    error('The sizes of the geometric matrices (mass, damping, stiffness) are not the same as added mass and damping matrices');
                end

                if (~hbcomp.isHB)
                    hbcomp.fbs = fbods;
                else
                    hbcomp.fbs = hydro;
                end
                
                % assume 1 device
                hbcomp.devCount = 1;
            end
        end
                
        function [iwavs] = get.IncWaves(hbcomp)
            % Incident waves
            iwavs = hbcomp.iwaves;
        end
        function [] = set.IncWaves(hbcomp, iwavs)
            if (~hbcomp.isHB)
                error('Can only set the incident waves for the computation when body is a HydroBody');
            end
            
            hbcomp.setIncWaves(iwavs);
        end
        
        function [wex] = get.Bodies(hbcomp)
            % Get the bodies in the array
            wex = hbcomp.fbs;
        end
        
        function [f] = get.Fex(hbcomp)
            % The hydrodynamic excitation force
            hbcomp.computeIfNot();
            
            f = hbcomp.fex;
        end
        
        function [f] = get.Ffk(hbcomp)
            % The Froude-Krylov force
            hbcomp.computeIfNot();
            
            f = hbcomp.ffk;
        end
        
        function [val] = get.DeviceCount(hbcomp)
            % The number of devices evaluated in the EnergyComp
            val = hbcomp.devCount;
        end
        function [] = set.DeviceCount(hbcomp, val)
            % The number of devices evaluated in the EnergyComp
            if ~isInt(val)
                error('The DeviceCount must be an integer');
            end
            hbcomp.devCount = val;
        end
        
        function [kfs, kfr] = KochinFuncs(hbcomp)
            
            if (~hbcomp.isHB)
                error('KochinFuncs can only computed when body is a HydroBody');
            end
            
            kfs(hbcomp.nT, hbcomp.nInc) = KochinFunc;
            kfr(hbcomp.nT, hbcomp.dof) = KochinFunc;
            
            dtm = hbcomp.fbs(1).DiffTM;
            
            for m = 1:hbcomp.nT
                M = hbcomp.fbs(1).Mlim(m);
                for n = 1:hbcomp.nInc
                    As = dtm{m}*hbcomp.iwaves(n).IncAmps(M, hbcomp.fbs(1).XYpos, m);
                    kfs(m, n) = KochinFunc('Coeffs', As.');
                end
            end
            
            Ar = hbcomp.fbs(1).RadCoefs;
            for m = 1:hbcomp.nT
                for n = 1:hbcomp.dof
                    kfr(m, n) = KochinFunc('Coeffs', Ar{m,n});
                end
            end 
        end
        
%         function [] = set.Bodies(hbcomp, fbs)
%             % Set the floating bodies in the array
% 
%             for n = 1:length(fbs)
%                 if (~isa(fbs(n), 'FloatingBody'))
%                     error('All floating bodies must be of type FloatingBody');
%                 end
%             end
%             
%             [m_, d_, k_] = IFreqDomComp.resizeMDK(fbs);
%             
%             [row col] = size(squeeze(hbcomp.a(1,:,:)));
%             
%             if (row ~= size(m_, 1) && col ~= size(m_, 2))
%                 error('The sizes of the geometric matrices (mass, damping, stiffness) are not the same as added mass and damping matrices');
%             end
%             
%             hbcomp.m = m_;
%             hbcomp.d = d_;
%             hbcomp.k = k_;
%             
%             hbcomp.fbs = fbs;
%             hbcomp.dof = IFreqDomComp.GetDoF(fbs);
%         end     

        function [wf, As] = WaveField(hbcomp, isarray, varargin)
            if (~hbcomp.isHB)
                error('Can only create a wavefield when body is a HydroBody');
            end
            
            hb = hbcomp.fbs(1);
            
            if (~isempty(hb.Rho))
                rho = hb.Rho;
            else
                rho = 1000;
            end
            
            DTM = hb.DiffTM;
            iwav = hbcomp.iwaves;
                                    
            for n = 1:hbcomp.nInc
                As = cell(hbcomp.nT, 1);                

                for m = 1:hbcomp.nT
                    Asm = DTM{m}*iwav(n).IncAmps(hb.Mlim(m), hb.XYpos, m);
                    As{m} = Asm.';
                end
                
                swaves = CirWaves('Out', hb.XYpos, As, hbcomp.t, hbcomp.h);
                swfs(n) = CirWaveField(rho, swaves, isarray, varargin{:});
                
                if (hbcomp.iwaves(n).IsPlane)
                    iwfs(n) = PlaneWaveField(rho, iwav(n), isarray, varargin{:});
                else
                    iwfs(n) = IncCirWaveField(rho, iwav(n), isarray, varargin{:});
                end
            end
            
            iwavefields = WaveFieldCollection(iwfs);
            swavefields = WaveFieldCollection(swfs);

            AR = hb.RadCoefs;
            
            for n = 1:hbcomp.dof
                Ar = AR(:,n);
                
                rwaves = CirWaves('Out', hb.XYpos, Ar, hbcomp.t, hbcomp.h);
                rwfs(n) = CirWaveField(rho, rwaves, isarray, varargin{:});
            end
            
            rwavefields = WaveFieldCollection(rwfs);
            
            wf = FBWaveField(iwavefields, swavefields, rwavefields);
            wf.BodyMotions = hbcomp.Motions;
        end
    end
      
    methods (Access = protected)
        function [] = computeIfNot(hbcomp) 
            if (hbcomp.isHB && ~hbcomp.isComp)
                fx = zeros(hbcomp.nT, hbcomp.nInc, hbcomp.DoF);
                hb = hbcomp.fbs(1);
                FTM = hb.FexTM;
                
                iwav = hbcomp.iwaves;
                for n = 1:hbcomp.nInc
                    for m = 1:hbcomp.nT
                        ftm = FTM{m};
                        M = hb.Mlim(m);
                        aI = iwav(n).IncAmps(M, hb.XYpos, m);
                        fx(m, n, :) = ftm*aI;
                    end
                end

                hbcomp.fex = fx;
                hbcomp.isComp = true;
            end
        end
        
        function [] = setAnyProp(obj, prop, val)
            obj.(prop) = val;
        end
    end
    
    methods (Static)
        function [newObj] = loadobj(a)
            if isstruct(a)
                newObj = feval(class(mfilename('class')));

                [currentUnset, loadedUnset] = CopyLoadHandleBase.trySetProps(newObj, a);
            else
                newObj = a;
            end
        end
    end
end