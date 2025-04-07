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
classdef WaveClimate < handle
    
    properties (Access = private)        
        specs;
        h;
        rho;
        name;
        hs;
        t02;
    end
    
    properties (Access = public)
        freqOcc;
    end
    
    properties (Dependent)
        WaveSpectra;
        H;
        T;
        Rho;
        Name;
    end
    
    methods 
        
        function [wc] = WaveClimate(spectra, freqOccur, varargin)
            if (nargin == 0)
                wc.specs = [];
                wc.freqOcc = [];
                wc.h = Inf;
                wc.rho = 1025;
            else
                if (isempty(varargin))
                    wc.h = Inf;
                    wc.rho = 1025;
                elseif (length(varargin) == 1)
                    wc.h = varargin{1};
                    wc.rho = 1025;
                else
                    wc.h = varargin{1};
                    wc.rho = varargin{2};
                end
                
                [opts, args] = checkOptions({{'Hs', 1}, {'T02', 1}}, varargin);
                if (opts(1))
                    wc.hs = args{1};
                else
                    wc.hs = [];
                end
                
                if (opts(2))
                    wc.t02 = args{2};
                else
                    wc.t02 = [];
                end
                
                WaveClimate.checkSpectra(spectra, freqOccur);
                
                wc.specs = spectra;
                wc.freqOcc = freqOccur;
            end
        end
                
        function [val] = get.Name(wc)
            val = wc.name;
        end
        function [] = set.Name(wc, val)
            wc.name = val;
        end
        
        function [ws] = get.WaveSpectra(wc)
            ws = wc.specs;
        end
        function [] = set.WaveSpectra(wc, ws)
            WaveClimate.checkSpectra(ws, wc.freqOcc)
            
            wc.specs = ws;
        end
        
        function [fo] = FreqOccurance(wc, varargin)
            opts = checkOptions({'Hours'}, varargin);
            ishrs = opts(1);
            if (ishrs)
                totHrsYr = 24*365;
                fo = totHrsYr*wc.freqOcc;
            else
                fo = wc.freqOcc;
            end
        end
        
        function [] = SetFreqOccurance(wc, fo, varargin)
            opts = checkOptions({'Hours'}, varargin);
            ishrs = opts(1);
            
            WaveClimate.checkSpectra(wc.specs, fo)
            if (ishrs)
                totHrsYr = 24*365;
                fo = fo./totHrsYr;
            end
            
            wc.freqOcc = fo;
        end
        
        function [rh] = get.Rho(wc)
            rh = wc.rho;
        end
        function [] = set.Rho(wc, rh)
            if (rh > 0)
                wc.rho = rh;
            else
                error('The density must be a positive number');
            end
        end 
        
        function [h_] = get.H(wc)
            % Get the water depeth
            h_ = wc.h;
        end
        function [] = set.H(wc, h_)
            % Set the water depth
            if (h_ > 0)
                wc.h = h_;
            else
                error('The depth must be a positive number');
            end
        end 
        
        function [t] = get.T(wc)
            t = [];
            if (~isempty(wc.specs))
                t = 1./wc.specs(1,1).Frequencies;
            end
        end
        
        function [hs] = Hs(wc, varargin)
            [opts] = checkOptions({'Intended'}, varargin);
            if (opts(1))
                hs = wc.hs;
            else
                [M, N] = size(wc.specs);
                hs = zeros(M, N);
                for m = 1:M
                    for n = 1:N
                        hs(m, n) = wc.specs(m,n).SigWaveHeight;
                    end
                end
            end
        end
        
        function [se] = Se(wc, varargin)
            [opts] = checkOptions({'Intended'}, varargin);
            if (opts(1))
                hs = wc.hs;
                t02 = wc.t02;
            else
                hs = wc.Hs('Intended');
                t02 = wc.T02('Intended');
            end
            
            H = hs*2*sqrt(2)/4; % Find wave height of that regular wave that has equivalent energy to the irregular wave of Hs.
            
            % Find Te from T02 values stored
            TtypeIn = 'T02'; TtypeOut = 'Te';
            Te = Bretschneider.ConverterT(t02, TtypeIn, TtypeOut);
            
            if size(hs,1) > 1 && size(hs,2) > 1
                L = repmat(1.5608*Te.*Te,size(hs,1),1); % Find wavelength of regular wave with period Te
                SeAll = H./L; % Find wave steepnesses
                Se = SeAll(:,1); % Pick out just one column - all should be identical if Hs matrix was generated using a single Se vector.
            else
                Harray = repmat(H',1,length(Te));
                Larray = repmat(1.5608*Te.*Te,length(hs),1);
                Se = Harray./Larray;
            end
            
            se = Se;
        end
        
        function [t02] = T02(wc, varargin)
            [opts] = checkOptions({'Intended'}, varargin);
            if (opts(1))
                t02 = wc.t02;
            else
                [M, N] = size(wc.specs);
                t02 = zeros(M, N);
                f = wc.specs(1,1).Frequencies;
                for m = 1:M
                    for n = 1:N
                        sp = wc.specs(m,n).Spectrum('Nondir');
                        m0 = abs(trapz(f,sp));
                        m2 = abs(trapz(f,f.^2.*sp));

                        t02(m, n) = sqrt(m0/m2);
                    end
                end
            end
        end
        
        function [M, N] = Size(wc)
            M = [];
            N = [];
            if (~isempty(wc.specs))
                [M, N] = size(wc.specs);
            end
        end
        
        function [Ef] = EnergyFlux(wc, varargin)
            [M, N] = size(wc.specs);
            Ef = cell(M, N);
            
            for m = 1:M
                for n = 1:N
                    ws = wc.specs(m, n);
                    Ef{m, n} = ws.EnergyFlux(wc.rho, wc. h, varargin{:});
                end
            end
        end
        
        function [spec] = AverageSpectrum(wc)
            [M, N] = size(wc.specs);
            f = wc.specs(1,1).Frequencies;
            dir = wc.specs(1,1).Directions;
%             if (isempty(dir))
%                 sp0 = zeros(1, length(f));
%             else
%                 sp0 = zeros(length(f), length(dir));
%             end
            
            sp0 = zeros(size(wc.specs(1,1).Spectrum));
            
            for m = 1:M
                for n = 1:N
                    sp = wc.specs(m,n).Spectrum;
                    fr = wc.freqOcc(m,n);
                    sp0 = sp0 + sp*fr;
                end
            end
            
%             if iscolumn(sp0)
%                 sp0 = sp0.';
%             end
            
            if (isempty(dir))
                spec = WaveSpectrum(sp0, f);
            else
                spec = WaveSpectrum(sp0, f, dir);
            end
        end
        
        function [ef] = AverageEnergyFlux(wc)
            spec = wc.AverageSpectrum;
            ef = spec.EnergyFlux(wc.rho, wc.h, 'total');
        end
        
        function [climI] = InterpolateTo(clim, Hs, T, varargin) % Interpolates clim into new Hs and T bins
            [opts, args] = checkOptions({{'typeSe',1}}, varargin);
            
            if opts(1)
                [t0M, hsM] = meshgrid(clim.t02, clim.Se);
                HsMatrix = Hs;
                Se = args{1};
                [TM, HsM] = meshgrid(T, Se);
            else
                [t0M, hsM] = meshgrid(clim.t02, clim.hs);
                [TM, HsM] = meshgrid(T, Hs);
            end
            
            freqOcc2 = interp2(t0M, hsM, clim.freqOcc, TM, HsM);
            
            freqOcc2(isnan(freqOcc2)) = 0;
            freqOcc2 = freqOcc2./sum(sum(freqOcc2));
            
            if isempty(clim.specs)
                climI = WaveClimate.MakeWaveClimate('Bretschneider', Hs, T,...
                    [], 'H', clim.h, 'Rho', clim.rho, 'T02');
            else
                climI = WaveClimate.MakeWaveClimate(class(clim.specs(1)), Hs, T,...
                    clim.specs(1).Frequencies, 'H', clim.h, 'Rho', clim.rho, 'T02');
            end
            climI.SetFreqOccurance(freqOcc2)
            climI.Name = clim.Name;
        end
        
        function [climo] = plus(clim, clim2)
            Hs1 = clim.Hs('Intended');
            T021 = clim.T02('Intended');
            Hs2 = clim2.Hs('Intended');
            T022 = clim2.T02('Intended');
            
            if (any(Hs1 ~= Hs2) || any (T021 ~= T022))
                error('To add wave climates, the matrix size must be the same');
            end
            
            freqOcc1 = clim.FreqOccurance;
            freqOcc2 = clim2.FreqOccurance;
            freqOcco = (freqOcc1 + freqOcc2)./2;
            
            climo = WaveClimate.MakeWaveClimate(class(clim.specs(1)), Hs1, T021,...
                clim.specs(1).Frequencies, 'Rho', clim.rho, 'T02');
            climo.SetFreqOccurance(freqOcco)
        end
        
        function [T0, Hs0] = PlotScatter(clim, varargin)
            [opts, args] = checkOptions({{'skip', 1}, {'hslim', 1}, ...
                {'percent'}, {'Tp'}, {'Te'}, {'EfWeight'}, {'occurance'}}, varargin);
            
            skip = 2;
            if opts(1)
                skip = args{1};
            end
            hslim = [];
            if opts(2)
                hslim = args{2};
            end
            percent = opts(3);
            useTp = opts(4);
            useTe = opts(5);
            efWeight = opts(6);
            occurance = opts(7);
            
            Hs = clim.Hs('Intended');
            T_ = clim.T02('Intended');
            if useTp
                T_ = T_./0.71;
            elseif useTe
                T_ = 1.206*T_;
            end
            
            if percent
                freqO = 100*clim.FreqOccurance;
                ylab = '% of year';
            elseif occurance
                freqO = 100*clim.FreqOccurance;
                ylab = 'Occurance (%)';
            else
                freqO = clim.FreqOccurance('hours');
                ylab = 'hours/year';
            end
            
            if efWeight
                Ef = cell2mat(clim.EnergyFlux('total'));
                if percent
                    freqO = Ef.*freqO./100./1000;
                    ylab = 'Power weighted (kW*frac)';
                else
                    freqO = Ef.*freqO./1000;
                    ylab = 'Power weighted (kWh)';
                end
            end
            
            iStopHs = length(Hs);
            if ~isempty(hslim)
                iStopHs = indexOf(Hs, hslim);
            end
            
            Hs = Hs(1:iStopHs);
            freqO = freqO(1:iStopHs, :);
            
            [~, inds] = maxnd(freqO);
            Hs0 = Hs(inds(1));
            T0 = T_(inds(2));
            
            indsHs = 1:skip:length(Hs);
            indsT = 1:skip:length(T_);
            
            plotScatter(T_, Hs, freqO, 'xinds', indsT, 'yinds', indsHs);

            if useTp
                xlabel('Tp (s)');
            elseif useTe
                xlabel('Te (s)');
            else
                xlabel('T02 (s)');
            end
            ylabel('Hs (m)');
            cb = colorbar;
            ylabel(cb, ylab);
        end
        
        function [AEP] = AnnualEnergyProduction(clim, Hsp, T02p, Pow)
            clim2 = clim.InterpolateTo(Hsp, T02p);
            freqH = clim2.FreqOccurance('hours');
            
            kWh = Pow.*freqH;
            
            AEP = sum(sum(kWh))/1000;
        end
    end
    
    methods (Static)
        function [wc] = MakeWaveClimate(type, Hs, T, f, varargin)
            [opts, args] = checkOptions({{'FreqOcc', 1}, {'H', 1}, {'Rho', 1}, {'T02'}, {'Tp'}, {'spread',4}, {'spread_logSpacing'}}, varargin);
            
            if size(Hs,1) > 1 && size(Hs,2) > 1
                typeSe = 1; % when wave steepness is used, there is a unique Hs for each Se-Te pairing
                nHs = size(Hs,1);
            else
                typeSe = 0;
                nHs = length(Hs);
            end
            nT = length(T);
            
            if (opts(1))
                freqOcc_ = args{1};
            else
                freqOcc_ = ones(nHs, nT);
            end
            
            if (opts(2))
                h_ = args{2};
            else
                h_ = Inf;
            end
            
            if (opts(3))
                rho_ = args{3};
            else
                rho_ = 1025;
            end
            
            if (opts(4))
                Ttype = 'T02';
            elseif (opts(5))
                Ttype = 'Tp';
            else
                Ttype = 'T02';
            end

            spread = opts(6);
            if spread
                s = args{6}{1}; % Spreading parameter for cos^(2s) distribution
                dirc = args{6}{2}; % Wave direction at centre of distribution
                Nb = args{6}{3}; % Set number of angles to use for wave components
                Nw = args{6}{4}; % Set number of frequencies to use for wave components

                logSpacing = opts(7); % Set whether to use log spacing or linear for frequencies

                % Catch errors
                if length(s) > 1
                    errorMsg = sprintf(['MoeCreateClimates is only set up to accommodate a single spreading parameter.\n'...
                        'The user should run a loop outside this function in order to create a 3D power matrix.']);
                    error(errorMsg)
                end
                if s < 1 || s > 100
                    warning('Spreading parameters outside range 1-100 are quite extreme. Are you sure this is correct?')
                end
                if length(dirc) > 1
                    error('User must specify a single value for the centre of the distribution.')
                end
                if length(Nb) > 1 || Nb <= 0 || rem(Nb,1) > 0
                    error('User must specify a single, positive, integer value for the centre of the distribution.')
                end
                if length(Nw) > 1 || Nw <= 0 || rem(Nw,1) > 0
                    error('User must specify a single, positive, integer value for the centre of the distribution.')
                end
            end

            if strcmpi(type, 'bretschneider')
                bs = true;
            elseif strcmpi(type, 'jonswap')
                bs = false;
            else
                error('Currently only Bretschneider and JONSWAP spectra are supported');
            end          
                        
            if isempty(f)
                specs_ = [];
            else
                specs_(nHs, nT) = Bretschneider;
                
                if spread
                    % With numWaveComponents = Nw*Nb, will want to trim
                    % spectrum first:
                    for m = 1:nHs
                        for n = 1:nT
                            if strcmp(Ttype,'T02')
                                T02_n = T(n);
                            end
                            if typeSe
                                Hs_m = Hs(m,n);
                            else
                                Hs_m = Hs(m);
                            end
                            testAngRange = [-pi/2:pi/200:pi/2]+dirc; testOmegaRange = 2*pi*f;
                            [~, AwTemp0] = WaveClimate.bretschneider(Hs_m, T02_n, testOmegaRange, 123, 'cos2s',s,dirc,testAngRange');

                            [wMaxVal wMaxInd] = max(max(abs(AwTemp0)));
                            angleInds = find(abs(AwTemp0(:,wMaxInd))>wMaxVal/100); % Find angle indices of above range that are worth considering in our spectrum.
                            betaForSpecAll(m,n,:) = linspace(0.99*testAngRange(min(angleInds)), 0.99*testAngRange(max(angleInds)), Nb);

                            [betaMaxVal betaMaxInd] = max(max(abs(AwTemp0),[],2));
                            freqInds = find(abs(AwTemp0(betaMaxInd,:))>betaMaxVal/100); % Find frequency indices of above range that are worth considering in our spectrum.
                            wForSpec = logspace(log10(1.01*testOmegaRange(min(freqInds))), log10(0.99*testOmegaRange(max(freqInds))), Nw);
                            fForSpecAll(m,n,:) = wForSpec/(2*pi);
                        end
                    end
                    minBeta = min(betaForSpecAll,[],'all'); % Find minimum beta from all individual trimmed spectra
                    maxBeta = max(betaForSpecAll,[],'all'); % Find maximum beta from all individual trimmed spectra
                    betaForSpec = linspace(minBeta, maxBeta, Nb);

                    minf = min(fForSpecAll,[],'all'); % Find minimum f from all individual trimmed spectra
                    maxf = max(fForSpecAll,[],'all'); % Find maximum f from all individual trimmed spectra
                    if logSpacing
                        fForSpec = logspace(log10(minf), log10(maxf), Nw);
                    else
                        fForSpec = linspace(minf, maxf, Nw);
                    end
                end

                for m = 1:nHs
                    for n = 1:nT
                        if bs
                            if spread
                                if typeSe
                                    specs_(m,n) = Bretschneider(Hs(m,n), T(n), 1./fForSpec, Ttype, varargin{:}, 'Cos2s', s, dirc, betaForSpec);
                                else
                                    specs_(m,n) = Bretschneider(Hs(m), T(n), 1./fForSpec, Ttype, varargin{:}, 'Cos2s', s, dirc, betaForSpec);
                                end
                            else
                                if typeSe
                                    specs_(m,n) = Bretschneider(Hs(m,n), T(n), 1./f, Ttype, varargin{:});
                                else
                                    specs_(m,n) = Bretschneider(Hs(m), T(n), 1./f, Ttype, varargin{:});
                                end
                            end
                        else
                            if typeSe
                                specs_(m,n) = JONSWAP(Hs(m,n), T(n), 1./f, Ttype, varargin{:});
                            else
                                specs_(m,n) = JONSWAP(Hs(m), T(n), 1./f, Ttype, varargin{:});
                            end
                        end
                    end
                end
            end
            
            if (strcmpi(Ttype, 'T02'))
                wc = WaveClimate(specs_, freqOcc_, h_, rho_, 'Hs', Hs, 'T02', T);
            else
                wc = WaveClimate(specs_, freqOcc_, h_, rho_, 'Hs', Hs);
            end
        end
        
        function [clim] = SumClimates(clims)
            Hs = clims(1).Hs('Intended');
            T02 = clims(1).T02('Intended');
            freqOcc = clims(1).FreqOccurance;
            
            Nclim = length(clims);
            
            for n = 2:Nclim
                Hsn = clims(n).Hs('Intended');
                T02n = clims(n).T02('Intended');
            
                if (any(Hs ~= Hsn) || any (T02 ~= T02n))
                    error('To add wave climates, the matrix size must be the same');
                end
            
                freqOccn = clims(n).FreqOccurance;
                freqOcc = freqOcc + freqOccn;
            end
            
            freqOcc = freqOcc./Nclim;

            clim = WaveClimate.MakeWaveClimate(class(clims(1).specs(1)), Hs, T02,...
                clims(1).specs(1).Frequencies, 'Rho', clims(1).rho, 'T02');
            clim.SetFreqOccurance(freqOcc)
        end
    end
     
    methods (Static, Access = private)
        function [] = checkSpectra(spectra, freqs)
            
            M = [];
            N = [];
            
            if (~isempty(spectra))
                [M, N] = size(spectra);
                f = spectra(1,1).Frequencies;
                beta = spectra(1,1).Directions;
                
                for m = 1:M
                    for n = 1:N
                        if (~isa(spectra(m,n), 'WaveSpectrum'));
                            error('All spectra must be of type WaveSpectrum');
                        end
                        
                        if (~all(spectra(m,n).Frequencies == f))
                            error('All spectra must have the same frequency range');
                        end
                        
                        if (~all(spectra(m,n).Directions == beta))
                            error('All spectra must have the same direction range');
                        end
                    end
                end
            end
            
            if (~isempty(freqs))
               [Mf, Nf] = size(freqs);
               if (~isempty(M))
                   if ((M ~= Mf) || (N ~= Nf))
                       error('The size of the array of spectra must be the same as the size of the array of frequencies of occurance');
                   end
               end
            end
        end

        function [S, A] = bretschneider(Hs, Tp, w, seed, varargin)

            [opts, args] = checkOptions({{'Cos2s', 3}}, varargin);
            spread = opts(1);

            wp = 2*pi./Tp;
            coef = 5/16*wp^4*Hs^2;
            Suni = coef./(w.^5).*exp(-5/4*(wp./w).^4);

            if spread
                s = args{1}{1};
                dirc = args{1}{2};
                dir = args{1}{3};

                G = cosSpectSpread(s, dirc, dir);
                S = Suni.*G;

                dbeta = dir(2:end) - dir(1:end-1);
                dbeta(end+1) = dbeta(end);
            else
                S = Suni;

                dbeta = 1;
            end

            if ~isempty(seed)
                rng(seed);
            end
            phase = (2*pi*rand(length(dbeta), length(w)));
            if iscolumn(w)
                phase = phase';
            end

            dw = w(2:end) - w(1:end-1);
            dw(end+1) = dw(end);

            A = (sqrt(2*S.*repmat(dw,length(dbeta),1).*repmat(dbeta,1,length(dw))).*exp(1i*phase));

        end

    end
    
end