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
classdef CreateClimates
    
    methods (Static)
        function [waveClim, Hs, T02] = EMEC(T, varargin)
            [waveClim, Hs, T02] = CreateClimates.readClimOriginLCOESS('emec', T, varargin{:});            
        end
        
        function [waveClim, Hs, T02] = BIMEP(T, varargin)
            [waveClim, Hs, T02] = CreateClimates.readClimOriginLCOESS('bimep', T, varargin{:});            
        end
        
        function [waveClim, Hs, T02] = YeuIsland(T, varargin)
            [waveClim, Hs, T02] = CreateClimates.readClimOriginLCOESS('yeu', T, varargin{:});            
        end
        
        function [waveClim, Hs, T02] = WaveHub(T, varargin)
            [waveClim, Hs, T02] = CreateClimates.readClimOriginLCOESS('hub', T, varargin{:});            
        end
        
        function [waveClim, Hs, T02] = DKNorthSea(T, varargin)
            [waveClim, Hs, T02] = CreateClimates.readClimOriginLCOESS('dk', T, varargin{:});            
        end
    end
    
    methods (Static, Access = private)
        function [waveClim, Hs, T02] = readClimOriginLCOESS(name, T, varargin)
            [opts, args] = checkOptions({{'HsLim', 1}}, varargin);
            hslim = inf;
            if opts(1)
                hslim = args{1};
            end
            
            switch name
                case 'emec'
                    clim = csvread([mwavePath '\Wave\WaveSpectrum\climates\emec.csv']);
                    climName = 'EMEC';
                case 'bimep'
                    clim = csvread([mwavePath '\Wave\WaveSpectrum\climates\bimep.csv']);
                    climName = 'BIMEP';
                case 'dk'
                    clim = csvread([mwavePath '\Wave\WaveSpectrum\climates\dk-northSea-pt2.csv']);
                    climName = 'North Sea (DK)';
                case 'yeu'
                    clim = csvread([mwavePath '\Wave\WaveSpectrum\climates\yeu-island.csv']);
                    climName = 'Yeu Island';
                case 'hub'
                    clim = csvread([mwavePath '\Wave\WaveSpectrum\climates\wave-hub.csv']);
                    climName = 'Wave Hub';
                otherwise
                    warning('Wave climate name not recongnized');
            end
            
            Hs = clim(2:end,1).';
            T02 = clim(1,2:end);
            freqOcc = clim(2:end, 2:end);
            
            [Hs, freqOcc] = CreateClimates.reduceMaxHs(Hs, freqOcc, hslim);

            waveClim = WaveClimate.MakeWaveClimate('bretschneider', Hs, T02, 1./T);
            waveClim.SetFreqOccurance(freqOcc, 'Hours')
            waveClim.Name = climName;
        end
        
        function [HsOut, freqOccOut] = reduceMaxHs(HsIn, freqOccIn, hslim)
            ihs = HsIn <= hslim;
            HsOut = HsIn(ihs);
            freqOccOut = freqOccIn(ihs,:);
        end
    end
end