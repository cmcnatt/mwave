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
            [opts, args] = checkOptions({{'HsLim', 1}}, varargin);
            if (opts(1))
                hslim = args{1};
            else
                hslim = Inf;
            end

            emec1 = csvread([mwavePath '\Wave\WaveSpectrum\climates\emec.csv']);
            
            Hs = emec1(2:end,1).';
            T02 = emec1(1,2:end);
            freqOcc = emec1(2:end, 2:end);
            
            ihs = Hs <= hslim;
            Hs = Hs(ihs);
            freqOcc = freqOcc(ihs,:);

            waveClim = WaveClimate.MakeWaveClimate('bretschneider', Hs, T02, 1./T);
            waveClim.SetFreqOccurance(freqOcc, 'Hours')
        end
    end
end