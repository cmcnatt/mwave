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

[opts, args] = checkOptions({{'SigFigCutoff', 1}, {'Movie', 1}, {'Mmax', 1}, ...
    {'AccTrim'}, {'Waves', 1}, {'Scattered', 1} {'Radiated', 1}, {'AxisSym',1}, ...
    {'XSym',1}, {'YSym',1}, {'XYSym', 1}, {'OldMeth'}}, varargin);
% 'Mmax' is an optional input that limits the number of coefficients to
% include in the diffraction transfer matrix and radiation coeffs
%
% 'SigFigCutoff', N tells the computation to only use the first N 
% significant figures. 'AccTrim' removes values on the ends that do not
% contribute significantly to reproducing the wave field. Without these
% options, when you compute arrays, you may get a warning message:
%
% Warning: Matrix is close to singular or badly scaled. Results may be 
% inaccurate. RCOND =  6.653070e-65. 
%
% This is because the diffraction transfer matrix has some small values,
% which are likely meaningless, and which makes the matrix solution hard to
% find. To get rid of these values, when you create the HydroBody, use, the
% options 'SigFigCutoff', 5, and 'AccTrim'. See example, Wamit_createHB_1

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

if (opts(8))
    symA = true;
    modes = args{8};
    motFuncs = modes.MotionFuncs;
else
    symA = false;
end

if (opts(9))
    symX = true;
    modes = args{9};
    motFuncs = modes.MotionFuncs;
else
    symX = false;
end

if (opts(10))
    symY = true;
    modes = args{10};
    motFuncs = modes.MotionFuncs;
else
    symY = false;
end

if (opts(11))
    symX = true;
    symY = true;
    modes = args{11};
    motFuncs = modes.MotionFuncs;
end

oldMeth = opts(12);

if (oldMeth && any([symA, symX, symY]))
    error('Old method not supported for symmetry');
end

if (symA && (symX || symY))
    error('Body cannot have both axis symmetry and symmetry in X or Y')
end

points = waveCP.FieldPoints;

h = waveCP.H;
beta = waveCP.IncWaveVals;
T = waveCP.T;
nT = length(T);

lam = round(IWaves.T2Lam(T,h));

MbA = length(beta);     % actual number of dicections

if (xor(symX, symY))
    Mb = 2*(MbA - 1);   % effective number of directions;
elseif (symX && symY)
    Mb = 4*(MbA - 1);
else
    Mb = MbA;
end

if (mod(Mb,2) == 1)
    Mb2 = (Mb-1)/2; 
else
    Mb2 = Mb/2-1;
end

if (Mmax < 0)
    if (symA)
        Mmax = 30;
    else
        Mmax = Mb2;
    end
else
    if (Mmax > Mb2)
        Mmax = Mb2;
    end
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
    % For each wave period (nT is the number of periods), find the
    % diffraction transfer matrix and the radiation amplitudes
    
    % 1) get scattering coeffs for the diffraction transfer matrix
    aSb = zeros(2*Mmax+1, MbA);   % aSb is an array of the scattered coeffs at each incident wave direction

    trimMs = zeros(MbA+dof,1);   % trimMs holds the cutoff M value

    for m = 1:MbA
        
        if (nT == 1)
            [r0, theta, z, etaq] = reshapeCirWFPoints(points, etaS{m});
        else
            [r0, theta, z, etaq] = reshapeCirWFPoints(points, etaS{nt, m});
        end
  
        [aS, eta0] = waveDecomp(k0(nt), r0, theta, z, h, etaq, L, Mmax);
        
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
        
        [aq, eta0] = waveDecomp(k0(nt), r0, theta, z, h, etaq, L, Mmax);
        
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

        [aq, trimMs(MbA+q)] = getTrimM(aq, eta0, sigFig, accTrim, theta, k0(nt), r0,...
            showMov, ['Radiated, \lambda/a = ' num2str(lam(nt)) ', Mode = ' num2str(q)], [movLoc 'rad_lam' num2str(lam(nt)) '_mod' num2str(q)]);

        AR{nt, q} = aq;
    end
    
    % Trim coefs
    M = max(trimMs);
    aSb = aSb(Mmax-M+1:Mmax+M+1,:);
    
    for q = 1:dof
        AR{nt, q} = AR{nt, q}(Mmax-M+1:Mmax+M+1);
    end
       
    % 3 - solve for each row of the diffraction transfer matrix B
    D = zeros(2*M+1, 2*M+1);

    if (oldMeth)
        % 2 - generate rotation matrix
        aIbeta = zeros(2*M+1, MbA);

        for m = -M:M
            for n = 1:MbA
                aIbeta(m+M+1, n) = exp(-1i*m*(pi/2 + beta(n)));
            end
        end
        aIbetaT = aIbeta.';
    
        for m = -M:M
            bm = aIbetaT\(aSb(m+M+1,:).');
            D(m+M+1,:) = bm.';
        end
    else
        expinpi = exp(1i*(-M:M)*pi/2);
        if (symA)
            for m = -M:M
                D(m+M+1,m+M+1) = expinpi(m+M+1)*aSb(m+M+1);
            end
        elseif (symX && ~symY)
            for m = -M:M
                for n = -M:M
                    Dmn = 1/(2*pi)*expinpi(n+M+1)*trapz(beta,(aSb(m+M+1,:).*exp(1i*n*beta)+exp(-1i*m*pi)*aSb(-m+M+1,:).*exp(-1i*n*beta)));
                    D(m+M+1, n+M+1) = Dmn;
                end
            end
        elseif (~symX && symY)
            for m = -M:M
                for n = -M:M
                    Dmn = 1/(2*pi)*expinpi(n+M+1)*trapz(beta,(aSb(m+M+1,:).*exp(1i*n*beta)+exp(1i*n*pi)*aSb(-m+M+1,:).*exp(-1i*n*beta)));
                    D(m+M+1, n+M+1) = Dmn;
                end
            end
        elseif (symX && symY)
            for m = -M:M
                for n = -M:M
                    Dmn = 1/(2*pi)*expinpi(n+M+1)*(cos((m+n)*pi) + 1)...
                        *trapz(beta,(aSb(m+M+1,:).*exp(1i*n*beta)+exp(-1i*m*pi)*aSb(-m+M+1,:).*exp(-1i*n*beta)));
                    D(m+M+1, n+M+1) = Dmn;
                end
            end
        else
            for m = -M:M
                X = ifft(aSb(m+M+1,:));

                bm = zeros(1, 2*M+1);
                bm(M+1) =  X(1);
                bm(M+2:2*M+1) = X(2:M+1);
                bm(M:-1:1) = X(MbA:-1:MbA-M+1);

                D(m+M+1,:) = expinpi.*bm;
            end
        end
    end
    
    if (sigFig > 0)
        D = roundSigFig(D, sigFig);
    end

    DTM{nt} = D;
    
    % Force Transfer Matrix
    if ((nT == 1) && (dof == 1))
        fbeta = Fex.';
    elseif (dof == 1)
        fbeta = squeeze(Fex(nt,:,:)).';
    else
        fbeta = squeeze(Fex(nt,:,:));
    end
    
    D = zeros(dof, 2*M+1);
    
    if (oldMeth)
        for q = 1:dof
            % force transfer matrix
            dq = aIbetaT\fbeta(:,q);
            D(q,:) = dq.';
        end
    else
        if (symA)
            for q = 1:dof
                switch motFuncs(q).ISym
                    case 1
                        D(q,M) = -1i/2*fbeta(q);        % m = -1
                        D(q,M+2) = 1i/2*fbeta(q);       % m = 1
                    case 2
                        warning('Symmtry forces assumes all 6 DOF are present for computation of sway and roll');
                        D(q,M) = -1/2*fbeta(q-1);       % m = -1
                        D(q,M+2) = -1/2*fbeta(q-1);      % m = 1
                    case 3
                        D(q,M+1) = fbeta(q);            % m = 0
                    case 4
                        warning('Symmtry forces assumes all 6 DOF are present for computation of sway and roll');
                        D(q,M) = 1/2*fbeta(q+1);        % m = -1
                        D(q,M+2) = 1/2*fbeta(q+1);      % m = 1
                    case 5
                        D(q,M) = -1i/2*fbeta(q);        % m = -1
                        D(q,M+2) = 1i/2*fbeta(q);       % m = 1
                    case 6
                end
            end
        elseif (symX && ~symY)
            for q = 1:dof
                if ((motFuncs(q).ISym == 1) || (motFuncs(q).ISym == 3) || (motFuncs(q).ISym == 5))
                    % surge, heave, pitch
                    for m = -M:M
                        Gqm = 1/pi*expinpi(m+M+1)*trapz(beta, fbeta(:,q).*cos(m*beta.'));
                        D(q, m+M+1) = Gqm;
                    end
                else
                    % sway, roll, yaw
                    for m = -M:M
                        Gqm = 1i/pi*expinpi(m+M+1)*trapz(beta, fbeta(:,q).*sin(m*beta.'));
                        D(q, m+M+1) = Gqm;
                    end
                end
            end
        elseif (~symX && symY)
            for q = 1:dof
                if ((motFuncs(q).ISym == 2) || (motFuncs(q).ISym == 3) || (motFuncs(q).ISym == 4))
                    % sway, heave, roll
                    evenMode = true;
                else
                    % surge, pitch, yaw
                    evenMode = false;
                end
                for m = -M:M
                    if (xor(evenMode, mod(m,2)))
                        Gqm = 1/pi*expinpi(m+M+1)*trapz(beta, fbeta(:,q).*cos(m*beta.'));
                        D(q, m+M+1) = Gqm;
                    else
                        Gqm = 1i/pi*expinpi(m+M+1)*trapz(beta, fbeta(:,q).*sin(m*beta.'));
                        D(q, m+M+1) = Gqm;
                    end
                end
            end
        elseif (symX && symY)
            for q = 1:dof
                if ((motFuncs(q).ISym == 1) || (motFuncs(q).ISym == 5))
                    % surge, pitch
                    for m = -M:M
                        if (mod(m,2))
                            Gqm = 2/pi*expinpi(m+M+1)*trapz(beta, fbeta(:,q).*cos(m*beta.'));
                            D(q, m+M+1) = Gqm;
                        end
                    end
                elseif ((motFuncs(q).ISym == 2) || (motFuncs(q).ISym == 4))
                    % sway, roll
                    for m = -M:M
                        if (mod(m,2))
                            Gqm = 2i/pi*expinpi(m+M+1)*trapz(beta, fbeta(:,q).*sin(m*beta.'));
                        	D(q, m+M+1) = Gqm;
                        end
                    end
                elseif (motFuncs(q).ISym == 3)
                    % heave
                    for m = -M:M
                        if (~mod(m,2))
                            Gqm = 2/pi*expinpi(m+M+1)*trapz(beta, fbeta(:,q).*cos(m*beta.'));
                            D(q, m+M+1) = Gqm;
                        end
                    end
                else
                    % yaw
                    for m = -M:M
                        if (~mod(m,2))
                            Gqm = 2i/pi*expinpi(m+M+1)*trapz(beta, fbeta(:,q).*sin(m*beta.'));
                        	D(q, m+M+1) = Gqm;
                        end
                    end
                end
            end
        else
            for q = 1:dof
                X = ifft(fbeta(:,q));

                dq = zeros(1, 2*M+1);
                dq(M+1) = X(1);
                dq(M+2:2*M+1) = X(2:M+1);
                dq(M:-1:1) = X(MbA:-1:MbA-M+1);
                D(q,:) = expinpi.*dq;
            end
        end
    end
    
    if (sigFig > 0)
        D = roundSigFig(D, sigFig);
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

function [Drnd] = roundSigFig(Din, sigFig)

Drnd = Din;
maxVal = max(max(abs(Din)));
expV = floor(log10(maxVal));
cutOffVal = 10^(expV-sigFig);
%izero = abs(Din) < cutOffVal;
Drnd = round(Drnd./cutOffVal)*cutOffVal;
%Drnd(izero) = 0;

end

