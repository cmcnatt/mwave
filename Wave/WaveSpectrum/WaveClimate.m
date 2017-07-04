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
        freqOcc;
        h;
        rho;
        name;
        hs;
        t02;
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
        
        function [climI] = InterpolateTo(clim, Hs, T)
            [t0M, hsM] = meshgrid(clim.t02, clim.hs);
            [TM, HsM] = meshgrid(T, Hs);
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
        
        function [] = PlotScatter(clim, varargin)
            [opts, args] = checkOptions({{'skip', 1}}, varargin);
            
            skip = 2;
            if opts(1)
                skip = args{1};
            end
            
            Hs = clim.Hs('Intended');
            T02 = clim.T02('Intended');
            freqO = clim.FreqOccurance('hours');
            
            indsHs = 1:skip:length(Hs);
            indsT = 1:skip:length(T02);
            
            plotScatter(T02, Hs, freqO, 'xinds', indsT, 'yinds', indsHs);

            xlabel('T02 (s)');
            ylabel('Hs (m)');
            cb = colorbar;
            ylabel(cb, 'hours/year');
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
            [opts, args] = checkOptions({{'FreqOcc', 1}, {'H', 1}, {'Rho', 1}, {'T02'}, {'Tp'}}, varargin);
            
            nHs = length(Hs);
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
            
            if strcmpi(type, 'bretschneider')
                bs = true;
            elseif strcmpi(type, 'jonswap')
                bs = false;
            else
                error('Currently only Bretschneider and JONSWAP specta are supported');
            end          
                        
            if isempty(f)
                specs_ = [];
            else
                specs_(nHs, nT) = Bretschneider;
                for m = 1:nHs
                    for n = 1:nT
                        if bs
                            specs_(m,n) = Bretschneider(Hs(m), T(n), 1./f, Ttype);
                        else
                            specs_(m,n) = JONSWAP(Hs(m), T(n), 1./f, Ttype);
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
    end
    
end