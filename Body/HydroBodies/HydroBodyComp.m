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
classdef HydroBodyComp < IHydroComp
    % Computes properties (motions, velocity, power) for a single or an 
    % array of floating bodies.
    
    properties (Access = private)
        fbs;
        isHB;
        fex;
    end
    
    properties (Dependent)
        IncWaves;
        Bodies;
        Fex;
    end
    
    methods
        
        function [hbcomp] = HydroBodyComp(hydro, varargin)
            % Constructor - takes either a single HydroBody of a HydroForces object and and a vector of
            % floating bodies
            
            hbcomp.isHB = false;
            
            if (isa(hydro, 'HydroBody'))
                hbcomp.isHB = true;
                if (isempty(varargin))
                    error('Second argument must be IWaves');
                elseif(~isa(varargin{1}, 'IWaves'));
                    error('Second argument must be IWaves');
                end
                iwavs = varargin{1};
            elseif (isa(hydro, 'HydroForces'))                
                if (isempty(varargin{1}))
                    error('Second input must be a FloatingBody');
                end
                fbods = varargin{1};
                if (length(varargin) > 1)
                    error('Can only set the incident waves for the computation when body is a HydroBody');
                else
                    nB = length(hydro.Beta);
                    wuns = ones(size(hydro.T));
                    for n = 1:nB
                        iwavs(n) = PlaneWaves(wuns, hydro.T, hydro.Beta(n)*wuns, hydro.H);
                    end
                end
            else
                error('Inputs must be either a HydroBody or a HydroForces object and and a vector of floating bodies');
            end
            
            if (~hbcomp.isHB)       
                for n = 1:length(fbods)
                    if (~isa(fbods(n), 'FloatingBody'))
                        error('Second input must be a FloatingBody')
                    end
                end
                
                hbcomp.initHydroParam(hydro.T, hydro.H, fbods);
                hbcomp.setIncWaves(iwavs);
                hbcomp.fex = hydro.Fex;
            else
                hbcomp.initHydroParam(hydro.T, hydro.H, hydro);
                hbcomp.setIncWaves(iwavs);
            end

            hbcomp.a = hydro.A;
            hbcomp.b = hydro.B;
            hbcomp.c = hydro.C;
            
            [row, col] = size(squeeze(hbcomp.a(1,:,:)));

            if (row ~= size(hbcomp.m, 1) && col ~= size(hbcomp.m, 2))
                error('The sizes of the geometric matrices (mass, damping, stiffness) are not the same as added mass and damping matrices');
            end

            if (~hbcomp.isHB)
                hbcomp.fbs = fbods;
            else
                hbcomp.fbs = hydro;
            end
        end
                
        function [iwavs] = get.IncWaves(hbcomp)
            % Incident waves
            iwavs = hbcomp.iwaves;
        end
        function [hbcomp] = set.IncWaves(hbcomp, iwavs)
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
%             [m_, d_, k_] = IHydroComp.resizeMDK(fbs);
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
%             hbcomp.dof = IHydroComp.GetDoF(fbs);
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
    end
end