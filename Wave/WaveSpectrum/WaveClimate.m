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
    end
    
    properties (Dependent)
        WaveSpectra;
        H;
        T;
        Hs;
        T02;
        Rho;
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
                WaveClimate.checkSpectra(spectra, freqOccur)

                wc.specs = spectra;
                wc.freqOcc = freqOccur;
            end
        end
                
        function [ws] = get.WaveSpectra(wc)
            ws = wc.specs;
        end
        function [wc] = set.WaveSpectra(wc, ws)
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
        function [wc] = set.Rho(wc, rh)
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
        function [wc] = set.H(wc, h_)
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
                t = w./wc.specs(1,1).Frequencies;
            end
        end
        
        function [hs] = get.Hs(wc)
            [M, N] = size(wc.specs);
            hs = zeros(M, N);
            for m = 1:M
                for n = 1:N
                    hs(m, n) = wc.specs(m,n).SigWaveHeight;
                end
            end
        end
        
        function [t02] = get.T02(wc)
            [M, N] = size(wc.specs);
            t02 = zeros(M, N);
            f = wc.specs(1,1);
            for m = 1:M
                for n = 1:N
                    sp = wc.specs(m,n).Spectrum('Nondir');
                    m0 = abs(trapz(f,sp));
                    m2 = abs(trapz(f,f.^2.*sp));

                    t02(m, n) = sqrt(m0/m2);
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
            
            if (~strcmpi(type, 'bretschneider'))
                error('Currently only bretschneider specta are supported');
            end          
                        
            specs_(nHs, nT) = Bretschneider;
            for m = 1:nHs
                for n = 1:nT
                    specs_(m,n) = Bretschneider(Hs(m), T(n), 1./f, Ttype);
                end
            end
            
            wc = WaveClimate(specs_, freqOcc_, h_, rho_);
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