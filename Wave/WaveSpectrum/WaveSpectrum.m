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
classdef WaveSpectrum <handle
    
    properties (Access = private)
        isdir;
        spectrum;
        frequencies;
        directions;
        cutoff;
        spreading;
    end
    
    properties (Dependent)
        IsDir;
        M0;
        SigWaveHeight;
        PeakPeriod;
        Spreading;
    end
    
    methods
        
        % Constructor
        function [spec] = WaveSpectrum(Spec, f, varargin)
            if (nargin > 0)                
                dir = spec.checkValues(Spec, f, varargin{:});
                                
                if (length(dir) == 0)
                    spec.isdir = 0;
                else
                    spec.isdir = 1;
                end

                spec.spectrum = Spec;
                spec.frequencies = f;
                spec.directions = dir;

                spec.cutoff = 0.01;
            end
        end
        
        % IsDirectional
        function [isDir] = get.IsDir(spec)
            isDir = spec.isdir;
        end
        
        % M0
        function [m0] = get.M0(spec)
            spec1 = spec.Spectrum('Nondir');

            del = spec.Deltas('Frequency');
            if isrow(spec1)
                spec1 = spec1.';
            end
            m0 = abs(sum(spec1.*del.'));
        end
        
        % SigWaveHeight
        function [Hs] = get.SigWaveHeight(spec)
            Hs = 4*sqrt(spec.M0);
        end

        % Spreading
        function [s] = get.Spreading(spec)
            s = spec.spreading;
        end

        function [] = set.Spreading(spec, s)
            spec.Spreading = s;
        end
        
        % PeakPeriod
        function [Tp] = get.PeakPeriod(spec)
            Spec = spec.Spectrum('Nondir');
            
            
            [mmm, imax] = max(Spec);
            fp = spec.frequencies(imax);
            Tp = 1./fp;
        end
        
        % Spectrum        
        function [Spec] = Spectrum(spec, varargin)
            [opts, args] = checkOptions({{'Nondir'}, {'WithEnergy'},...
                {'Period'}, {'Wavelength', 1}}, varargin);
            
            nondir = opts(1);
            withEnergy = opts(2);
            period = opts(3);
            wavelen = opts(4);
            if wavelen
                h = args{4};
            end
            
            if period || wavelen
                
                ind = 1;
                args = {};
                if nondir
                    args{ind} = 'Nondir';
                    ind = ind + 1;
                end
                
                if withEnergy
                    args{ind} = 'WithEnergy';
                    ind = ind + 1;
                end
                    
                a = spec.Amplitudes(args{:});
                f = spec.Frequencies(args{:});
                T = 1./f;
                if (period)
                    del = computeDelta(T);
                elseif (wavelen)
                    lam = IWaves.T2Lam(T, h);
                    del = computeDelta(lam);
                end
                
                if spec.isdir
                    Spec = spec.getNonDensitySpec;
                else
                    df = spec.Deltas('Frequency');
                    Spec = spec.spectrum.*df;
                end
                
                Spec = abs(Spec./del);
            else
                if nondir
                    if (spec.isdir)
                        del = spec.Deltas('Direction');
                        %del = del.'*ones(1, length(spec.frequencies));
                        Spec = (spec.spectrum*del.');
                    else
                        Spec = spec.spectrum;
                    end

                    if (withEnergy)
                        [~, ifreq] = spec.findFreqWithEnergy(spec.cutoff);
                        Spec = Spec(ifreq);
                    end
                else
                    Spec = spec.spectrum;
                    if (withEnergy)
                        [~, ifreq] = spec.findFreqWithEnergy(spec.cutoff);

                        if (spec.isdir)
                            [~, idir] = spec.findDirWithEnergy(spec.cutoff);
                            Spec = Spec(ifreq, idir);
                        else
                            Spec = Spec(ifreq);
                        end
                    end
                end
            end
        end
        
        % Frequencies
        function [freq] = Frequencies(spec, varargin)
            
            withEnergy = 0;
            
            if (~isempty(varargin))
                if (strcmp(varargin{1}, 'WithEnergy'))
                    withEnergy = 1;
                end
            end
                        
            if (withEnergy)
                freq = spec.findFreqWithEnergy(spec.cutoff);
            else
                freq = spec.frequencies;
            end
        end
        
        % Directions
        function [dir] = Directions(spec, varargin)
            if (spec.isdir)
                withEnergy = 0;

                if (~isempty(varargin))
                    if (strcmp(varargin{1}, 'WithEnergy'))
                        withEnergy = 1;
                    end
                end

                if (withEnergy)
                    dir = spec.findDirWithEnergy(spec.cutoff);
                else
                    dir = spec.directions;
                end
            else
                dir = [];
            end
        end
        
        % Deltas
        function [del1, del2] = Deltas(spec, type)
            switch (type)
                case 'Frequency'
                    del1 = computeDelta(spec.frequencies);
                case 'Direction'
                    if (length(spec.directions) > 1)
                        del1 = computeDelta(spec.directions);
                    else
                        del1 = 1;
                    end
                case 'Both'
                    del1 = computeDelta(spec.frequencies);
                    if (length(spec.directions) > 1)
                        del2 = computeDelta(spec.directions);
                    else
                        del2 = 1;
                    end
                otherwise
                    error('Deltas takes an argument of ''Frequency'', ''Direction'' of ''Both''');
            end
        end
        
        % Amplitudes
        function [a] = Amplitudes(spec, varargin)
            [opts, args] = checkOptions({{'Nondir'}, {'WithEnergy'}, ...
                {'RandPhase'}, {'RandSeed', 1}}, varargin);
            
            nondir = opts(1);
            withEnergy = opts(2);
            randPhase = opts(3);
            seed = [];
            if opts(4)
                randPhase = true;
                seed = args{4};
            end
                              
            if (nondir || ~spec.isdir)
                Spec = spec.Spectrum('Nondir');
                df = spec.Deltas('Frequency');
                a = sqrt(2*Spec.*df);
                
                if (withEnergy)
                    [~, ifreq] = spec.findFreqWithEnergy(spec.cutoff);

                    a = a(ifreq);
                end
            else
                a = sqrt(2*spec.getNonDensitySpec);

                if (withEnergy)
                    [~, ifreq] = spec.findFreqWithEnergy(spec.cutoff);
                    [~, idir] = spec.findDirWithEnergy(spec.cutoff);
                    
                    a = a(ifreq, idir);
                end
            end
            
            if (imag(a(1,1)) ~= 0)
                a = imag(a);
            end
            
            if (randPhase)
                if ~isempty(seed)
                    rng(seed);
                end
                a = a.*exp(1i*2*pi*rand(size(a)));
            end
        end

        function [waves] = GetPlaneWaves(spec, h, varargin)
            f = spec.Frequencies(varargin{:});

            beta = spec.Directions(varargin{:});
            if (isempty(beta))
                beta = 0;
            end
            nb = length(beta);
            
            a = spec.Amplitudes(varargin{:});

            if (nb > 1)
                for o = 1:length(beta)
                    waves(o) = PlaneWaves(a(:,o), 1./f, beta(o), h);
                end
            else
                waves = PlaneWaves(a, 1./f, 0, h);
            end
        end
        
        function [Ef] = EnergyFlux(spec, rho, h, varargin)
            opts = checkOptions({{'Total'}, {'Density'}}, varargin);
            tot = opts(1);
            density = opts(2);
            
            nf = length(spec.Frequencies(varargin{:}));
            if (density)
                df = spec.specs(1,1).Deltas('Frequency');
                if (df(1) < 0)
                    df = -df;
                end
            end

            iwaves = spec.GetPlaneWaves(h, varargin{:});
            nb = length(iwaves);

            ef = zeros(nf, nb);
            if (nb > 1)
                for o = 1:nb
                    efo = iwaves(o).EnergyFlux(rho);
                    ef(:,o) = efo;
                end
            else
                ef = iwaves.EnergyFlux(rho);
            end

            if (tot)
                Ef = sum(sum(ef));
            else
                if(density)
                    Ef = ef./df;
                else
                    Ef = ef;
                end
            end
        end
        
        function [H, T] = EnergyWave(spec, rho, h)
            Ef = spec.EnergyFlux(rho, h, 'Total');
            if ~isa(spec, 'Bretschneider')
                error('EnergyWave computation not set up for non-Bretschnider wave');
            end
            T = Bretschneider.ConverterT(spec.PeakPeriod, 'Tp', 'Te');
            wav = PlaneWaves(1, T, 0, h);

            a = sqrt(2*Ef/(rho*IWaves.G*wav.Cg));
            H = 2*a;
        end
        
        function [] = Redistribute(spec, type, vals, varargin)
            
            S = spec.spectrum;
            
            % Creating a new interpolation surface, assumes that the directions wrap around 
            dir = spec.directions;
            
            nMid = floor(length(dir)/2);
            maxDir = dir(1) + dir(end);
            dir1 = dir(1:nMid);
            dir2 = dir(nMid+1:end);
            S1 = S(:,1:nMid);
            S2 = S(:,nMid+1:end);
            
            dirb = [dir2-maxDir;dir;dir1+maxDir]; 
            Sb = [S2,S,S1];
            
            f = spec.frequencies;
            
            [Dir, F] = meshgrid(dirb,f);
            
            switch (type)
                case 'Frequency'
                    fn = vals;
                    dirn = spec.directions;                    
                case 'Direction'
                    fn = spec.frequencies;
                    dirn = vals;   
                case 'Both'
                    fn = vals;
                    dirn = varargin{1};
                otherwise
                    error('Redistribute takes an argument of ''Frequency'', ''Direction'' of ''Both''');
            end
            
            [Dirn, Fn] = meshgrid(dirn, fn);
            
            Sn = interp2(Dir, F, Sb, Dirn, Fn);
       
            spec.frequencies = fn;
            spec.directions = dirn;
            spec.spectrum = Sn;
        end
        
        function [Spec2] = Copy(spec, varargin)
            if (~isempty(varargin))
                if (strcmp(varargin{1}, 'WithEnergy'))
                    Spec2 = WaveSpectrum(spec.Spectrum('WithEnergy'), spec.Frequencies('WithEnergy'), spec.Directions('WithEnergy'));
                end
            else
                Spec2 = WaveSpectrum(spec.spectrum, spec.frequencies, spec.directions);
            end
        end
        
    end
    
    methods (Access = protected)
        function [] = setSpectrum(spec, specin, f, varargin)
            dir = spec.checkValues(specin, f, varargin{:});
            
            if (length(dir) == 0)
                spec.isdir = 0;
            else
                spec.isdir = 1;
            end
            
            spec.spectrum = specin;
            spec.frequencies = f;
            spec.directions = dir;
            if length(varargin) == 2
                spec.spreading = varargin{2};
            end
            
            spec.cutoff = 0.01;
        end
    end
    
    methods (Access = private)
        
        function [dir] = checkValues(spec, specin, f, varargin)
            if (~isempty(varargin) > 0)
                dir = varargin{1};
            else
                dir = [];
            end
                
            if (isvector(f))
                if (isrow(f))
                    f = f.';
                end
            else
                error('The frequncies input must be a vector');
            end
                                
            nf = length(f);
            ndir = length(dir);
                
            if (ndir == 0)
                spec.isdir = 0;

                if (isvector(specin))
                    if (length(specin) ~= nf)
                        error('The input spectrum must be the same size as the frequencies');
                    end
                    if (isrow(specin))
                        specin = specin.';
                    end
                else
                    error('For the nondirectional case, the spectrum input must be a vector');
                end
            else
                spec.isdir = 1;

                if (isvector(dir))
                    if (iscolumn(dir))
                        dir = dir.';
                    end
                else
                    error('The directions input must be a vector');
                end

                [row, col] = size(specin);
                if (ndir ~= col || nf ~= row)
                    error('The spectrum must be a matirx of size Nfreq x Ndir');
                end
            end
        end
        
        function [freq, ifreq] = findFreqWithEnergy(spec, cutoff)
             maxVal = max(max(spec.spectrum));
             if (spec.isdir)
                 ifreq = max(spec.spectrum, [], 2) >= cutoff*maxVal;
             else
                 ifreq = spec.spectrum >= cutoff*maxVal;
             end
             freq = spec.frequencies(ifreq);
        end
        
        function [dir, idir] = findDirWithEnergy(spec, cutoff)
             if (spec.isdir)
                 maxVal = max(max(spec.spectrum));
                 idir = max(spec.spectrum) >= cutoff*maxVal;
                 dir = spec.directions(idir);
             else
                 dir = [];
                 idir = [];
             end
        end
        
        function [Spec] = getNonDensitySpec(spec)
            [df, dd] = spec.Deltas('Both');
                
            [Dd, Df] = meshgrid(dd, df);
            
            Spec = spec.spectrum.*Df.*Dd;
        end
    end
    
end