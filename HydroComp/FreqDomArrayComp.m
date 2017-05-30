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
classdef FreqDomArrayComp < IFreqDomComp
    % Computes array interactions - Fex, A, B for an array of floating bodies 
    % using Kagemoto and Yue's intecation theory. Floating bodies must be of 
    % type HydroBody - that is, they have a diffraction transfer matrix, a 
    % force transfer matrix, and their own added mass and damping.
    %
    % Properties (motions, velocity, power) are computed with the parent 
    % FreqDomComp class.
    
    properties (Access = private)
        omega;
        k0;
        hbs;
        Nhb;
        isBTmatComp;
        compTime;
        Dmat;
        Tmat;
        Mmax;
        Mval;
        fex;
        As;
        Ar;
    end
    
    properties (Dependent)
        IncWaves;
        Bodies;
        BodXY;
        BodAng;
        CompTime;
        Fex;
    end
    
     methods
        
         function [hacomp] = FreqDomArrayComp(hydroBs, iwaves)            
            % Constructor - takes a vector of HydroBody objects
            
            if (~isa(hydroBs(1), 'HydroBody'))
                error('Inputs must be a vector of HydroBody objects');
            end
            
            hacomp.Nhb = length(hydroBs);
            
            for n = 2:hacomp.Nhb
                if (hydroBs(n) ~= hydroBs(1))
                    error('All HydroBodies must have the same wave periods, water depth, and water density');
                end
            end
            
            hacomp.initHydroParam(hydroBs(1).T, hydroBs(1).H, hydroBs, false, 0);
            hacomp.omega = 2*pi./hacomp.T;
            hacomp.k0 = IWaves.SolveForK(hacomp.omega, hacomp.h);
            hacomp.hbs = hydroBs;
            
            hacomp.setIncWaves(iwaves)
            
            hacomp.computeMmax();
            
            hacomp.isComp = false;
            hacomp.isBTmatComp = zeros(hacomp.nT,1);
            hacomp.compTime = 0;
            hacomp.Dmat = cell(hacomp.nT,1);
            hacomp.Tmat = cell(hacomp.nT,1);
            hacomp.Mval = zeros(hacomp.nT,1);
         end
         
         function [iwavs] = get.IncWaves(hbcomp)
            % Incident waves
            iwavs = hbcomp.iwaves;
        end
        function [hbcomp] = set.IncWaves(hbcomp, iwavs)
            hbcomp.setIncWaves(iwavs);
        end
         
         function [wex] = get.Bodies(hacomp)
            % Get the bodies in the array
            wex = hacomp.hbs;
         end  
         
         function [xy] = get.BodXY(hacomp)
             % Get the XY position of the bodies in the array
             xy = zeros(hacomp.Nhb,2);
             for n = 1:hacomp.Nhb
                 xy(n,:) = hacomp.hbs(n).XYpos;
             end
         end
         function [] = set.BodXY(hacomp, xy)
             % Set the XY position of the bodies in the array
             [M, N] = size(xy);
             if ((M ~= hacomp.Nhb) || (N ~= 2))
                 error('The xy position input must be an array of size (Number of Bodies)x2');
             end
             
             for n = 1:hacomp.Nhb
                 hacomp.hbs(n).XYpos = xy(n,:);
             end        
             
             hacomp.isComp = false;
             hacomp.isBTmatComp = zeros(hacomp.nT,1);
             hacomp.compTime = 0;
         end
         
         function [ang] = get.BodAng(hacomp)
             % Get the rotation angles of the bodies in the array
             ang = zeros(hacomp.Nhb, 1);
             for n = 1:hacomp.Nhb
                 ang(n) = hacomp.hbs(n).Angle;
             end
         end
         function [] = set.BodAng(hacomp, ang)
             % Set the XY position of the bodies in the array
             [M, N] = size(ang);
             if ((M ~= hacomp.Nhb) || (N ~= 1))
                 error('The angle input must be an array of size (Number of Bodies)x2');
             end
             
             for n = 1:hacomp.Nhb
                 if ((ang(n) > 360) || (ang(n) < 0))
                     error('The angle input must be in degrees and must be between 0 and 360');
                 end
                 hacomp.hbs(n).Ang = ang(n);
             end        
             
             hacomp.isComp = false;
             hacomp.isBTmatComp = zeros(hacomp.nT,1);
             hacomp.compTime = 0;
         end
         
         function [ti] = get.CompTime(hacomp)
             ti = hacomp.compTime;
         end
         
         function [f] = get.Fex(hacomp)
            % The hydrodynamic excitation force
            hacomp.computeIfNot();
            
            f = hacomp.fex;
        end
    
         function [wf] = WaveField(hacomp, isarray, varargin)
             % Get the wave field
             hacomp.computeIfNot();
             hb = hacomp.hbs(1);
            
             if (~isempty(hb.Rho))
                 rho = hb.Rho;
             else
                 rho = 1000;
             end
             
             if (~isempty(hb.H))
                 h = hb.H;
             else
                 error('To compute wave fields, water depth cannot be empty');
             end

             zwf = ZeroWaveField(rho, hacomp.t, h, isarray, varargin{:});
             
             % hacomp.As = cell(hacomp.nT, hacomp.nInc, hacomp.Nhb);
             
             for n = 1:hacomp.nInc
                 swfsn = zwf;
                 for m = 1:hacomp.Nhb
                     Asnm = squeeze(hacomp.As(:,n,m));
                     for l = 1:hacomp.nT
                         Asnm{l} = Asnm{l}.';
                     end
                     swaves = CirWaves('Out', hacomp.hbs(m).XYpos, Asnm, hacomp.t, hacomp.h);
                     swfsn = swfsn + CirWaveField(rho, swaves, isarray, varargin{:});
                 end
                 swfs(n) = swfsn;
                 if (hacomp.iwaves(n).IsPlane)
                     iwfs(n) = PlaneWaveField(rho, hacomp.iwaves(n), isarray, varargin{:});
                 else
                     iwfs(n) = IncCirWaveField(rho, hacomp.iwaves(n), isarray, varargin{:});
                 end
            end
            
            iwavefields = WaveFieldCollection(iwfs);
            swavefields = WaveFieldCollection(swfs);
            
            rwfs(hacomp.dof, 1) = WaveField;
                        
            for n = 1:hacomp.dof 
                rwfsn = zwf;
                
                for m = 1:hacomp.Nhb
                    Arnm = squeeze(hacomp.Ar(:, n, m));
                    for l = 1:hacomp.nT
                        Arnm{l} = Arnm{l}.';
                    end
                
                    rwaves = CirWaves('Out', hacomp.hbs(m).XYpos, Arnm, hacomp.t, hacomp.h);
                    rwfsn = rwfsn + CirWaveField(rho, rwaves, isarray, varargin{:});
                end
                rwfs(n) = rwfsn;
            end
            
            rwavefields = WaveFieldCollection(rwfs);
            
            wf = FBWaveField(iwavefields, swavefields, rwavefields);
            wf.BodyMotions = hacomp.Motions;
         end
     end
     
     methods (Access = protected)
         
        function [] = computeIfNot(hacomp)
            if (~hacomp.isComp)
                tic;
                hacomp.computeArrayHydro();
                hacomp.compTime = toc;
            end
        end
        
    end
    
     methods (Static)

     end
     
     methods (Access = private)
        
         function [] = computeArrayHydro(hacomp)
             
             hacomp.checkBodySpacing();
             
             hacomp.computeDiffraction();
             
             hacomp.computeRadiation();
             
             hacomp.computeHydrostatic();
             
             hacomp.isComp = true;
         end
         
         function [] = computeMmax(hacomp)
             
             hacomp.Mmax = zeros(hacomp.nT, 1);
             
             for m = 1:hacomp.nT
                 mmax = 0;

                 for n = 1:hacomp.Nhb
                     M = hacomp.hbs(n).Mlim(m);
                     if (M > mmax)
                         mmax = M;
                     end
                 end

                 hacomp.Mmax(m) = mmax;
             end
         end
         
         function [] = checkBodySpacing(hacomp)
             % check spacing requirement - center of body i must not be
             % inside the circumscribed circle of body j

             bstart = 1;
             
             while (bstart ~= hacomp.Nhb)
                 for n = bstart:hacomp.Nhb
                     bi = hacomp.hbs(n);
                     Ri = bi.Rcir;
                     for m = bstart+1:hacomp.Nhb
                         if (m ~= n)
                             bj = hacomp.hbs(m);
                             Rj = bj.Rcir;
                             Lij = sqrt((bi.XYpos(1) - bj.XYpos(1))^2 + (bi.XYpos(2) - bj.XYpos(2))^2);

                             if ((Lij < Ri ) || (Lij < Rj))
                                 error('center of body i must not be inside the circumscribed circle of body j');
                             end
                         end
                     end
                 end
                 bstart = bstart + 1;
             end
         end
         
         function [] = computeDiffraction(hacomp)
             % compute the Fex for each wave period and direction, compute
             % As for each body for each wave period and direction
             
             hacomp.fex = zeros(hacomp.nT, hacomp.nInc, hacomp.dof);
             hacomp.As = cell(hacomp.nT, hacomp.nInc, hacomp.Nhb);
                          
             for m = 1:hacomp.nT
                 %M = hacomp.Mmax(m);
                 for n = 1:hacomp.nInc
                     
                     [aI, aSmn, M] = hacomp.incScatAmps(m, n);
                     hacomp.As(m, n, :) = aSmn;
                     
                     dofStart = 1;
                     for l = 1:hacomp.Nhb
                         dofStop = dofStart + hacomp.hbs(l).DoF - 1;
                         
                         Ml = hacomp.hbs(l).Mlim(m);
                         G = hacomp.hbs(l).FexTM{m};
                         
                         aIl = aI{l};
                         
                         if (Ml < M)
                             aIl = aIl(M-Ml+1:M+Ml+1);
                         end
                         
                         if (Ml > M)
                             G = G(:, Ml-M+1:Ml+M+1);
                         end
                         
                         f = G*aIl;
                         
                         hacomp.fex(m, n, dofStart:dofStop) = f;
                         
                         dofStart = dofStop + 1;
                     end
                 end
             end
         end
         
         function [] = computeRadiation(hacomp)
             % compute the A and B
             hacomp.a = zeros(hacomp.nT, hacomp.dof, hacomp.dof);
             hacomp.b = zeros(hacomp.nT, hacomp.dof, hacomp.dof);
             hacomp.Ar = cell(hacomp.nT, hacomp.dof, hacomp.Nhb);
             
             for l = 1:hacomp.nT
             
                 %M = hacomp.Mmax(l);
                 colDof = 1;
                 ldof1 = 1;
                 omeg = hacomp.omega(l);
                 Z = zeros(hacomp.dof, hacomp.dof);
                 
                 for m = 1:hacomp.Nhb
                     
                     dofm = hacomp.hbs(m).DoF;
                     ldof2 = ldof1 + dofm - 1;
                     
                     Am = squeeze(hacomp.hbs(m).A(l,:,:));
                     Bm = squeeze(hacomp.hbs(m).B(l,:,:));
                     Zm = 1i*omeg*Am + Bm;
                     
                     for q = 1:dofm
                         [aIlmq, aRlmq, M] = hacomp.incRadAmps(l, m, q);
                         
                         zq = zeros(hacomp.dof, 1);
                         zq(ldof1:ldof2) = Zm(:,q);
                         
                         qdof1 = 1;
                         
                         for n = 1:hacomp.Nhb
                             
                             Mn = hacomp.hbs(n).Mlim(l);
                             G = hacomp.hbs(n).FexTM{l};
                             aIlmqn = aIlmq{n};
                             
                             if (Mn < M)
                                 aIlmqn = aIlmqn(M-Mn+1:M+Mn+1);
                             end
                             
                             if (Mn > M)
                                 G = G(:, Mn-M+1:Mn+M+1);
                             end
                         
                             dofn = hacomp.hbs(n).DoF;
                             qdof2 = qdof1 + dofn - 1;
                             
                             zn = 1i/omeg*G*aIlmqn;         % This is -1/(i*omega)*G*aI
                             zq(qdof1:qdof2) = zq(qdof1:qdof2) + zn;
                             
                             qdof1 = qdof2 + 1;
                         end

                         Z(:, colDof) = zq;
                         
                         hacomp.Ar(l, colDof, :) = aRlmq;
                         
                         colDof = colDof + 1;
                     end
                     
                     ldof1 = ldof2 + 1;
                 end
                 hacomp.a(l,:,:) = imag(Z)./omeg;
                 hacomp.b(l,:,:) = real(Z);
             end
         end
         
         function [] = computeHydrostatic(hacomp)
             
             hacomp.c = zeros(hacomp.dof, hacomp.dof);
             dof1 = 1;
             
             for n = 1:hacomp.Nhb
                 dof = hacomp.hbs(n).DoF;
                 dof2 = dof1 + dof - 1;
                 
                 hacomp.c(dof1:dof2, dof1:dof2) = hacomp.hbs(n).C;
                 
                 dof1 = dof2 + 1;
             end
         end
         
         function [Ai, As, M] = incScatAmps(hacomp, iT, iwav)            
    
             [D, T, M] = hacomp.getDTMmats(iT);
             %M = hacomp.Mmax(iT);
             Nelem = hacomp.Nhb*(2*M + 1);

             aIa = zeros(Nelem, 1);
             
             incWav = hacomp.iwaves(iwav);
             
             for n = 1:hacomp.Nhb
                 % make incident wave amp matrix
                 n1 = (n-1)*(2*M+1)+1;
                 n2 = n1 + 2*M;
                 
                 aIa(n1:n2) = incWav.IncAmps(M, hacomp.hbs(n).XYpos, iT);
             end

             [aS, aI] = hacomp.solveScatAmps(Nelem, D, T, aIa);
             
             As = cell(hacomp.Nhb, 1);
             Ai = cell(hacomp.Nhb, 1);
             
             for n = 1:hacomp.Nhb
                 n1 = (n-1)*(2*M+1)+1;
                 n2 = n1 + 2*M;
                 
                 As{n} = aS(n1:n2);
                 Ai{n} = aI(n1:n2);
             end
         end
         
         function [Ai, Ar, M] = incRadAmps(hacomp, iT, iBod, iMode)
             
             [D, T, M] = hacomp.getDTMmats(iT);
             %M = hacomp.Mmax(iT);
             Nelem = hacomp.Nhb*(2*M + 1);

             aIa = zeros(Nelem, 1);
             
             for n = 1:hacomp.Nhb
                 n1 = (n-1)*(2*M+1)+1;
                 n2 = n1 + 2*M;
                 
                 % make radiated incident wave amp matrix
                 if (n ~= iBod)
                     aIa(n1:n2) = FreqDomArrayComp.rAmps(iT, M, hacomp.hbs(n), hacomp.hbs(iBod), iMode, hacomp.k0(iT));
                 end
             end
                         
             [aS, aI] = hacomp.solveScatAmps(Nelem, D, T, aIa);
                          
             aRi = hacomp.hbs(iBod).RadCoefs{iT, iMode};
             
             Mi = length(aRi);
             Mi = (Mi - 1)/2;
             
             aRiM = zeros(2*M+1, 1);
             
             if (M > Mi)
                 aRiM(M+1-Mi:M+1+Mi) = aRi;
             else
                 aRiM(:) = aRi(Mi+1-M:Mi+1+M);
             end
             
             Ar = cell(hacomp.Nhb, 1);
             Ai = cell(hacomp.Nhb, 1);
             
             for n = 1:hacomp.Nhb
                 n1 = (n-1)*(2*M+1)+1;
                 n2 = n1 + 2*M;
                 
                 if (n == iBod)
                     Ar{n} = aS(n1:n2) + aRiM;
                 else
                     Ar{n} = aS(n1:n2);
                 end
                 Ai{n} = aI(n1:n2);
             end
         end
         
         function [aS, aI] = solveScatAmps(hacomp, Nelem, D, T, aIa)
             A = eye(Nelem, Nelem) - D*T;
             b = D*aIa;
             
             aS = A\b;
             aI = aIa + T*aS;
         end
         
         function [D, T, M] = getDTMmats(hacomp, iT)
             
             if (~hacomp.isBTmatComp(iT))
                 
                 Mm = hacomp.Mmax(iT);
                 M = Mm;
                 
                 Nelem = hacomp.Nhb*(2*M + 1);

                 D = zeros(Nelem, Nelem);
                 T = zeros(Nelem, Nelem);

                 for m = 1:hacomp.Nhb
                     
                     m1 = (m-1)*(2*Mm+1)+1;
                     mCent = m1 + Mm;
                     m2 = m1 + 2*Mm;

                     % make large DFT - D
                     Dm = hacomp.hbs(m).DiffTM{iT};
                     
                     [Mn, ~] = size(Dm);
                     Mn = (Mn - 1)/2;

                     mB1 = mCent - Mn;
                     mB2 = mCent + Mn;

                     D(mB1:mB2, mB1:mB2) = Dm;

                     for n = 1:hacomp.Nhb
                         if (n ~= m)
                             n1 = (n-1)*(2*Mm+1) + 1;
                             n2 = n1 + 2*Mm;

                             Tmn = CirWaves.BasisTH(Mm, hacomp.hbs(m).XYpos, hacomp.hbs(n).XYpos, hacomp.k0(iT));

                             T(m1:m2, n1:n2) = Tmn.';
                         end
                     end
                 end

                 hacomp.Dmat{iT} = D;
                 hacomp.Tmat{iT} = T;
                 hacomp.Mval(iT) = M;
                 hacomp.isBTmatComp(iT) = true;
             else
                 D = hacomp.Dmat{iT};
                 T = hacomp.Tmat{iT};
                 M = hacomp.Mval(iT);
             end
         end
     end
        
     
     methods (Static)
                  
         function [aIj] = rAmps(iT, M, incBod, radBod, radMode, k0)
             
             Tij = CirWaves.BasisTH(M, incBod.XYpos, radBod.XYpos, k0);
             
             air = radBod.RadCoefs{iT, radMode};
             
             Mj = length(air);
             Mj = (Mj - 1)/2;
             
             airM = zeros(2*M+1, 1);
             if (M > Mj)
                 airM(M+1-Mj:M+1+Mj) = air;
             else
                 airM(:) = air(Mj+1-M:Mj+1+M);
             end

             aIj = (Tij.')*airM;
         end
     end
end