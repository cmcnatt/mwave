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
classdef IWavesUT < matlab.unittest.TestCase
    % Unit test for IWaves Class

    methods (Test)
       
        function testPlane1(testCase)
            
            a = 1;
            t = 1;
            beta = pi/2;
            h = 5;
            
            wc = PlaneWaves(a, t, beta, h);
            
            testCase.verifyEqual(wc.A, a);
            testCase.verifyEqual(wc.T, t);
            testCase.verifyEqual(wc.Beta, beta);
            testCase.verifyEqual(wc.H, h);
            testCase.verifyEqual(wc.Epsilon, 0);
            
            epsilon = pi;
            
            wc = PlaneWaves(a, t, beta, h, epsilon);
            testCase.verifyEqual(wc.Epsilon, epsilon);
        end
        
        function testPlane2(testCase)
            
            a = [1 0.5];
            t = [1 2];
            beta = [pi/2 pi/2];
            h = 5;
            
            wc = PlaneWaves(a, t, beta, h);
            
            testCase.verifyEqual(wc.A, a.');
            testCase.verifyEqual(wc.T, t.');
            testCase.verifyEqual(wc.Beta, beta.');
            testCase.verifyEqual(wc.H, h);
            testCase.verifyEqual(wc.Epsilon, [0; 0]);
            testCase.verifyEqual(wc.UniformDir, true);
           
            
            a = [1 0.5];
            t = [1 2];
            beta = [pi/2 pi/4];
            h = 5;
            
            wc = PlaneWaves(a, t, beta, h);
            
            testCase.verifyEqual(wc.A, a.');
            testCase.verifyEqual(wc.T, t.');
            testCase.verifyEqual(wc.Beta, beta.');
            testCase.verifyEqual(wc.H, h);
            testCase.verifyEqual(wc.Epsilon, [0; 0]);
            testCase.verifyEqual(wc.UniformDir, false);
        end
        
        function testProperties(testCase)
            a = 1;
            t = 1;
            beta = pi/2;
            h = 5;
            
            wc = PlaneWaves(a, t, beta, h);
            
            f = 1./t;
            testCase.verifyEqual(wc.F, f);
            
            omega = 2*pi*f;
            testCase.verifyEqual(wc.Omega, omega);
            
            k = IWaves.SolveForK(omega, h);
            testCase.verifyEqual(wc.K, k);
            
            lambda = 2*pi/k;
            testCase.verifyEqual(wc.Lambda, lambda);
            
            c = lambda/t;
            testCase.verifyEqual(wc.C, c);
            
            testCase.verifyEqual(wc.DepthType, 'Deep');
            
            cg = 0.5*c;
            testCase.verifyEqual(wc.Cg, cg, 'AbsTol', 1e-10);
            
            rho = 1000;
            g = IWaves.G;
            e = 0.5*rho*g*a^2;
            testCase.verifyEqual(wc.Energy(rho), e, 'AbsTol', 1e-10);
            
            ef = e*cg;
            testCase.verifyEqual(wc.EnergyFlux(rho), ef, 'AbsTol', 1e-10);
            
            % set properties
            a2 = 0.5;
            wc.A = a2;
            testCase.verifyEqual(wc.A, a2);
            try
                wc.A = [1 2];
                % If there is NO error, test fails
                testCase.verifyFail();
            catch
                % If there is an error, test passes
            end
            try
                wc.A = -1;
                % If there is NO error, test fails
                testCase.verifyFail();
            catch
                % If there is an error, test passes
            end
            
            epsilon = pi/4;
            wc.Epsilon = epsilon;
            testCase.verifyEqual(wc.Epsilon, epsilon);
            try
                wc.Epsilon = [1 2];
                % If there is NO error, test fails
                testCase.verifyFail();
            catch
                % If there is an error, test passes
            end
            
            beta = pi/4;
            wc.Beta = beta;
            testCase.verifyEqual(wc.Beta, beta);
            try
                wc.Beta = [1 2];
                % If there is NO error, test fails
                testCase.verifyFail();
            catch
                % If there is an error, test passes
            end
            
            t = 0.4;
            wc.T = t;
            testCase.verifyEqual(wc.T, t);
            try
                wc.T = [1 2];
                % If there is NO error, test fails
                testCase.verifyFail();
            catch
                % If there is an error, test passes
            end
            
            f = 0.4;
            wc.F = f;
            testCase.verifyEqual(wc.F, f);
            testCase.verifyEqual(wc.T, 1./f);
            testCase.verifyEqual(wc.Omega, 2*pi*f);
            k = IWaves.SolveForK(2*pi*f, wc.H);
            testCase.verifyEqual(wc.K, k);
            testCase.verifyEqual(wc.Lambda, 2*pi/k);
            testCase.verifyEqual(wc.C, 2*pi*f/k, 'AbsTol', 1e-10)
            try
                wc.F = [1 2];
                % If there is NO error, test fails
                testCase.verifyFail();
            catch
                % If there is an error, test passes
            end
            
            omega = pi;
            wc.Omega = omega;
            testCase.verifyEqual(wc.Omega, omega);
            testCase.verifyEqual(wc.T, 2*pi/omega);
            testCase.verifyEqual(wc.F, omega/(2*pi));
            k = IWaves.SolveForK(omega, wc.H);
            testCase.verifyEqual(wc.K, k);
            testCase.verifyEqual(wc.Lambda, 2*pi/k);
            testCase.verifyEqual(wc.C, omega/k, 'AbsTol', 1e-10);
            try
                wc.Omega = [1 2];
                % If there is NO error, test fails
                testCase.verifyFail();
            catch
                % If there is an error, test passes
            end
            
            k = 0.5;
            wc.K = k;
            testCase.verifyEqual(wc.K, k);
            testCase.verifyEqual(wc.Lambda, 2*pi/k);
            omega = IWaves.SolveForOmega(k, wc.H);
            testCase.verifyEqual(wc.Omega, omega);
            testCase.verifyEqual(wc.T, 2*pi/omega);
            testCase.verifyEqual(wc.F, omega/(2*pi));
            testCase.verifyEqual(wc.C, omega/k, 'AbsTol', 1e-10);
            try
                wc.K = [1 2];
                % If there is NO error, test fails
                testCase.verifyFail();
            catch
                % If there is an error, test passes
            end
            
            lambda = 3;
            wc.Lambda = lambda;
            testCase.verifyEqual(wc.Lambda, lambda);
            testCase.verifyEqual(wc.K, 2*pi/lambda);
            omega = IWaves.SolveForOmega(2*pi/lambda, wc.H);
            testCase.verifyEqual(wc.Omega, omega);
            testCase.verifyEqual(wc.T, 2*pi/omega, 'AbsTol', 1e-10);
            testCase.verifyEqual(wc.F, omega/(2*pi));
            testCase.verifyEqual(wc.C, omega*lambda/(2*pi), 'AbsTol', 1e-10);
            try
                wc.Lmabda = [1 2];
                % If there is NO error, test fails
                testCase.verifyFail();
            catch
                % If there is an error, test passes
            end
            
        end
        
        function testStaticMethods(testCase)
            
            t = 5;
            h = 20;
            rho = 1000;
            
            lam = 38.8976;
            
            testCase.verifyEqual(IWaves.T2Lam(t, h), lam, 'AbsTol', 1e-4);
            testCase.verifyEqual(IWaves.Lam2T(lam, h), t, 'AbsTol', 1e-4);
            
            k = 2*pi/lam;
            testCase.verifyEqual(IWaves.SolveForK(2*pi/t, h), k, 'AbsTol', 1e-4);
            testCase.verifyEqual(IWaves.SolveForOmega(k, h), 2*pi/t, 'AbsTol', 1e-4);
            
            ef = 1.9457e+04;
            testCase.verifyEqual(IWaves.UnitEnergyFlux(rho, t, h), ef, 'AbsTol', 1);
            
            ks = [0.1615   0.1081   0.2887   0.4542   0.6155];
            testCase.verifyEqual(IWaves.SolveForK(2*pi/t, h, 'Evanescent', 4), ks, 'AbsTol', 1e-4);
            
            h = Inf;
            k = (2*pi/t)^2./IWaves.G;
            testCase.verifyEqual(IWaves.SolveForK(2*pi/t, h), k, 'AbsTol', 1e-4);
            testCase.verifyEqual(IWaves.SolveForOmega(k, h), 2*pi/t, 'AbsTol', 1e-4);
            
            h = 0;
            try
                IWaves.SolveForK(2*pi/t, h);
                % If there is NO error, test fails
                testCase.verifyFail();
            catch
                % If there is an error, test passes
            end
            
            try
                IWaves.SolveForOmega(k, h);
                % If there is NO error, test fails
                testCase.verifyFail();
            catch
                % If there is an error, test passes
            end
            
        end
        
        function cirWavesTrans(testCase)
            % Scattered waves from j incident to points in i's coordinate system
            % Only accurate inside the circle around Ci of radius Lij
            
            %           |       Cj
            %           | P
            %           |   Ci
            %           |
            %---------------------------
            %           |   
            %           |       
            %           |
            %           |
            
            P = [1 3];
            
            x = -10:0.2:10;
            y = -10:0.2:10;
            Np = length(x);
            
            [X, Y] = meshgrid(x, y);
            
            k0 = 2*pi/2;
            
            Cj = [5 5];
            Pj = P - Cj;
            rPj = sqrt(Pj(1)^2 + Pj(2)^2);
            thetaPj = atan2(Pj(2), Pj(1));
            Xj = X - Cj(1);
            Yj = Y - Cj(2);
            Rj = sqrt(Xj.^2 + Yj.^2);
            Thetaj = atan2(Yj, Xj);
            
            Ci = [2 2];
            Pi = P - Ci;
            rPi = sqrt(Pi(1)^2 + Pi(2)^2);
            thetaPi = atan2(Pi(2), Pi(1));
            Xi = X - Ci(1);
            Yi = Y - Ci(2);
            Ri = sqrt(Xi.^2 + Yi.^2);
            Thetai = atan2(Yi, Xi);
            
            Mj = 3;
            Psij = zeros(2*Mj+1, 1);
            psij = zeros(2*Mj+1, Np, Np);            
            for m = -Mj:Mj
                Psij(m+Mj+1) = besselh(m, 2, k0*rPj)*exp(1i*m*thetaPj);
                psij(m+Mj+1, :, :) = besselh(m, 2, k0*Rj).*exp(1i*m*Thetaj);
            end
            
            Mi = 20;
            Psii = zeros(2*Mi+1, 1);
            psii = zeros(2*Mi+1, Np, Np);
            for m = -Mi:Mi
                Psii(m+Mi+1) = besselj(m, k0*rPi).*exp(1i*m*thetaPi);
                psii(m+Mi+1, :, :) = besselj(m, k0*Ri).*exp(1i*m*Thetai);
            end
            
            Tij = CirWaves.BasisTH(Mi, Ci, Cj, k0);
            
            Psij2 = Tij*Psii;
            Psij2 = Psij2(Mi-Mj+1:Mi+Mj+1);
            psij2 = zeros(2*Mj + 1, Np, Np);
            
            for m = 1:Np
                for n = 1:Np
                    psij21 = Tij*squeeze(psii(:, m, n));
                    psij2(:, m, n) = psij21(Mi-Mj+1:Mi+Mj+1);
                end
            end
            
            L = sqrt((Ci(1) - Cj(1))^2 + (Ci(2) - Cj(2))^2);
            cirLix = Ci(1) + L*cos(0:pi/20:2*pi);
            cirLiy = Ci(2) + L*sin(0:pi/20:2*pi);

            figure;
            
            testCase.verifyEqual(Psij(Mj+1), Psij2(Mj+1), 'AbsTol', 1e-6);
            
            subplot(Mj+1, 2, 1);
            pcolor(X, Y, real(squeeze(psij(Mj+1,:,:))));
            set(gca, 'clim', [-1 1], 'xtick', [], 'ytick', []);
            shading flat; axis equal; axis tight;
            hold on;
            plot([Cj(1) Ci(1)], [Cj(2) Ci(2)],'w', 'MarkerSize',8, 'Marker','*','LineStyle','none');
            plot(cirLix, cirLiy, 'w');
            ylabel('m = 0');
            title('Re(\Psi_j)');
            
            subplot(Mj+1, 2, 2);
            pcolor(X, Y, real(squeeze(psij2(Mj+1,:,:))));
            set(gca, 'clim', [-1 1], 'xtick', [], 'ytick', []);
            shading flat; axis equal; axis tight;
            hold on;
            plot([Cj(1) Ci(1)], [Cj(2) Ci(2)],'w', 'MarkerSize',8, 'Marker','*','LineStyle','none');
            plot(cirLix, cirLiy, 'w');
            title('Re(Tij\Psi_i)');
            
            for m = 1:Mj
                mp = Mj+1+m;
                mn = Mj+1-m;
                testCase.verifyEqual(Psij(mp), Psij2(mp), 'AbsTol', 1e-6);
                testCase.verifyEqual(Psij(mn), Psij2(mn), 'AbsTol', 1e-6);
                
                subplot(Mj+1,2,2*m+1);
                %psijPlot = 0.5*real(squeeze(psij(mp,:,:)) + (-1)^m*squeeze(psij(mn,:,:)));
                psijPlot = real(squeeze(psij(mp,:,:)));
                pcolor(X, Y, psijPlot);
                set(gca, 'clim', [-1 1], 'xtick', [], 'ytick', []);
                shading flat; axis equal; axis tight;
                hold on;
                plot([Cj(1) Ci(1)], [Cj(2) Ci(2)],'w', 'MarkerSize',8, 'Marker','*','LineStyle','none');
                plot(cirLix, cirLiy, 'w');
                %ylabel(['m = |' num2str(m) '|']);
                ylabel(['m = ' num2str(m)]);

                subplot(Mj+1,2,2*m+2);
                %psij2Plot = 0.5*real(squeeze(psij2(mp,:,:)) + (-1)^m*squeeze(psij2(mn,:,:)));
                psij2Plot = real(squeeze(psij2(mp,:,:)));
                pcolor(X, Y, psij2Plot);
                set(gca, 'clim', [-1 1], 'xtick', [], 'ytick', []);
                shading flat; axis equal; axis tight;
                hold on;
                plot([Cj(1) Ci(1)], [Cj(2) Ci(2)],'w', 'MarkerSize',8, 'Marker','*','LineStyle','none');
                plot(cirLix, cirLiy, 'w');
            end

            
        end
    end
    
end

