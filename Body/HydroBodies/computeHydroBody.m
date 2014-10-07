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
function [hydroB, eta0aSw, eta0aRw] = computeHydroBody(waveCP, hydroF, floatB, varargin)
% TODO: document this man

[opts, args] = checkOptions({{'SigFigCutoff', 1}, {'Movie', 1}, {'Mmax', 1}, {'AccTrim'}, {'Waves', 1}, {'Scattered', 1} {'Radiated', 1}}, varargin);

sigFig = args{1};
movLoc = args{2};
Mmax = args{3};

if (opts(2))
    showMov = 'Movie';
else
    showMov = [];
    movLoc = [];
end

if (opts(4))
    accTrim = 'AccTrim';
else
    accTrim = [];
end

retWaves = opts(5);
if (retWaves)
    iTw = args{5};
    if (opts(6))
        mbSw = args{6};
    else
        mbSw = -1;
    end
    
    if (opts(7))
        qRw = args{7};
    else
        qRw = -1;
    end
    
    eta0aSw = cell(length(mbSw),5);
    iSw = 1;
    eta0aRw = cell(length(qRw),5);
    iRw = 1;
end

points = waveCP.FieldPoints;

h = waveCP.H;
Beta = waveCP.IncWaveVals;
T = waveCP.T;
nT = length(T);

lam = round(IWaves.T2Lam(T,h));

Mb = length(Beta);
if (mod(Mb,2) == 1)
    Mb2 = (Mb-1)/2;
else
    Mb2 = Mb/2-1;
end

if ((Mmax > 0) && (Mmax < Mb2))
    Mb2 = Mmax;
end

L = 0;

k0 = IWaves.SolveForK(2*pi./T, h);

% Diffraction transfer matrix
etaS = waveCP.Elevation('Scattered');
DTM = cell(nT, 1);

% Force transfer matrix
dof = hydroF.DoF;
Fex = hydroF.Fex;
FTM = cell(nT, 1);

% Radiation coeffs
etaR = waveCP.Elevation('Radiated');
AR = cell(nT, dof);

for nt = 1:nT
    % Get scattering coeffs
    aSb = zeros(2*Mb2+1, Mb);

    trimMs = zeros(Mb+dof,1);

    for m = 1:Mb
        if (nT == 1)
            [r0, theta, z, etaq] = reshapeCirWFPoints(points, etaS{m});
        else
            [r0, theta, z, etaq] = reshapeCirWFPoints(points, etaS{nt, m});
        end
  
        [aS, eta0] = waveDecomp(k0(nt), r0, theta, z, h, etaq, L, Mb2);
        
        if (m == 1)
            showMov1 = showMov;
        else
            showMov1 = {};
        end
        
        if (retWaves)
            if (nt == iTw)
                if (any(m == mbSw))
                    eta0aSw{iSw, 1} = eta0;
                    eta0aSw{iSw, 2} = aS;
                    eta0aSw{iSw, 3} = theta;
                    eta0aSw{iSw, 4} = k0(nt);
                    eta0aSw{iSw, 5} = r0;
                    iSw = iSw+1;
                end
            end
        end
       
        [aS, trimMs(m)] = getTrimM(aS, eta0, sigFig, accTrim, theta, k0(nt), r0, ...
            showMov1, ['Scatter, \lambda/a = ' num2str(lam(nt))], [movLoc 'scat_lam' num2str(lam(nt))]);

        aSb(:,m) = aS.';
    end
    
    % Get radiation coeffs
    for q = 1:dof
        % radiation coeffs
        if (nT == 1)        
            [r0, theta, z, etaq] = reshapeCirWFPoints(points, etaR{q});
        else
            [r0, theta, z, etaq] = reshapeCirWFPoints(points, etaR{nt, q});
        end
        
        [aq, eta0] = waveDecomp(k0(nt), r0, theta, z, h, etaq, L, Mb2);
        
        if (retWaves)
            if (nt == iTw)
                if (any(q == qRw))
                    eta0aRw{iRw, 1} = eta0;
                    eta0aRw{iRw, 2} = aq;
                    eta0aRw{iRw, 3} = theta;
                    eta0aRw{iRw, 4} = k0(nt);
                    eta0aRw{iRw, 5} = r0;
                    iRw = iRw+1;
                end
            end
        end

        [aq, trimMs(Mb+q)] = getTrimM(aq, eta0, sigFig, accTrim, theta, k0(nt), r0,...
            showMov, ['Radiated, \lambda/a = ' num2str(lam(nt)) ', Mode = ' num2str(q)], [movLoc 'rad_lam' num2str(lam(nt)) '_mod' num2str(q)]);

        AR{nt, q} = aq;
    end
    
    % Trim coefs
    M = max(trimMs);
    aSb = aSb(Mb2-M+1:Mb2+M+1,:);
    
    for q = 1:dof
        AR{nt, q} = AR{nt, q}(Mb2-M+1:Mb2+M+1);
    end
    

    % 2 - generate rotation matrix
    aIbeta = zeros(2*M+1, Mb);

    for m = -M:M
        for n = 1:Mb
            aIbeta(m+M+1, n) = exp(-1i*m*(pi/2 + Beta(n)));
        end
    end
    
    aIbetaT = aIbeta.';

    % 3 - solve for each row of the diffraction transfer matrix B
    B = zeros(2*M+1, 2*M+1);

    for m = -M:M
        bm = aIbetaT\(aSb(m+M+1,:).');
        B(m+M+1,:) = bm.';
    end

    DTM{nt} = B;
    
    % Force Transfer Matrix
    if ((nT == 1) && (dof == 1))
        fbeta = Fex.';
    elseif (dof == 1)
        fbeta = squeeze(Fex(nt,:,:)).';
    else
        fbeta = squeeze(Fex(nt,:,:));
    end
    
    D = zeros(dof, 2*M+1);

    for q = 1:dof
        % force transfer matrix
        dq = aIbetaT\fbeta(:,q);
        D(q,:) = dq.';
    end
    
    FTM{nt} = D;
end
    
hydroB = HydroBody(floatB, hydroF, DTM, FTM, AR);
end

function [Atrim, trimM] = getTrimM(A, eta0, sigFigs, varargin)

[opts, args] = checkOptions({{'AccTrim', 3}, {'Movie', 2}}, varargin);

accTrim = opts(1);
if (accTrim)
    theta = args{1}{1};
    k0 = args{1}{2};
    r0 = args{1}{3};
end

showMov = opts(2);
if (showMov)
    movLoc = args{2}{2};
    name = args{2}{1};
end

M = length(A);
M = (M-1)/2;

trimM = M;

if (sigFigs > 0)
    maxVal = max(abs(eta0));
    expV = floor(log10(maxVal));
    cutOffVal = 10^(expV-sigFigs);

    for m = M:-1:1
        Am = abs(A(M+1+m));
        Anm = abs(A(M+1-m));

        if (all(cutOffVal >= Am) && all(cutOffVal >= Anm))
            trimM = m-1;
        else
            break;
        end
    end
    Atrim = zeros(1, 2*M+1);
    Atrim(M+1-trimM:M+1+trimM) = A(M+1-trimM:M+1+trimM);
    Atrim = round(Atrim/cutOffVal)*cutOffVal;
else
    Atrim = A;
end

if (accTrim)    
    H0 = besselh(0, 2, k0*r0);
    Hm = besselh(1:M, 2, k0*r0);
    Hnm = (-1).^(1:M).*Hm;

    H = [Hnm(M:-1:1), H0, Hm];

    if (showMov)
        figure;
    end

    for Mt = 0:trimM
        eta = zeros(size(theta));

        for m = -Mt:Mt
            eta = eta + A(M+1+m)*H(M+1+m)*exp(1i*m*theta);
        end

        err = abs(eta - eta0)./abs(eta0);

        if (showMov)
            plot(theta, abs(eta0), 'k');
            hold on;
            plot(theta, abs(eta));

            title({name, ['M = ' num2str(Mt) ', Mean error = ' num2str(mean(err))]});

            mov(Mt+1) = getframe(gcf);
            cla;
        end

        if (mean(err) < 0.01)
            trimM = Mt;
            break;
        end
    end
    
    if (showMov)
        close;
        movie2avi(mov, movLoc, 'compression','None', 'fps', 2);
    end
end
end

