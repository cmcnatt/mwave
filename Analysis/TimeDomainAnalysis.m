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
classdef TimeDomainAnalysis < handle
    % Holds and processes time-domain results, which could come from an
    % experiment or a simulation
    
    properties (Access = protected)    
        nSig;
        motDof;
        ptoDof;
        motions;
        ptoKinematic;
        ptoDynamic;
        waveSigs;
        waveAmps;
        fmotions;
        fptoKinematic;
        fptoDynamic;
        fwaveSigs;
        smotions;
        sptoKinematic;
        sptoDynamic;
        swaveSigs;
        motTime;
        ptoTime;
        waveTime;
        motFreq;
        ptoFreq;
        waveFreq;
        motTimeLims;
        ptoTimeLims;
        waveTimeLims;
        wgPos;
        waveBeta;
        meanPos;
    end
    
    properties (Dependent)
        MotionDoF;
        PtoDoF;
        NSignals;
        MeanBodyPosition;
    end
    
    methods
        
        function [tdan] = TimeDomainAnalysis()
        end
        
        function [val] = get.MotionDoF(tda)
            val = tda.motDof;
        end
        
        function [val] = get.PtoDoF(tda)
            val = tda.ptoDof;
        end
        
        function [val] = get.NSignals(tda)
            val = tda.nSig;
        end
        
        function [] = set.MeanBodyPosition(tda, val)
            tda.meanPos = val;
        end
        function [val] = get.MeanBodyPosition(tda)
            val = tda.meanPos;
        end
        
        function [] = SetMotions(tda, dofs, time, motions)
            tda.setValues('motions', dofs, time, motions);
        end
        
        function [] = SetPtoKinematic(tda, dofs, time, ptoKin)
            tda.setValues('ptoKinematic', dofs, time, ptoKin);
        end
        
        function [] = SetPtoDynamic(tda, dofs, time, ptoKin)
            tda.setValues('ptoDynamic', dofs, time, ptoKin);
        end
                
        function [] = SetWaves(tda, wgPos, time, waves, varargin)
            [Nwg, ~] = size(wgPos);
            
            isAmp = false;
            if isempty(time)
                isAmp = true;
                if isempty(varargin)
                    beta = 0;
                else
                    beta = varargin{1};
                end
            end
            
            if Nwg == 1 && ndims(waves) == 2
                [Nsig, Nta] = size(waves);
                wavesTemp = waves;
                waves = zeros(Nsig, 1, Nta);
                waves(:, 1, :) = wavesTemp;
            elseif Nwg > 1 && ndims(waves) == 2
               [NwgSig, Nta] = size(waves);
               wavesTemp = waves;
               waves = zeros(1, NwgSig, Nta);
               waves(1, :, :) = wavesTemp; 
            end
            
            [Nsig, NwgSig, Nta] = size(waves);
            
            err = 0;
            if isempty(tda.nSig)
                tda.nSig = Nsig;
            elseif (tda.nSig ~= Nsig)
                err = 1;
            end
            
            if NwgSig ~= Nwg
                err = 1;
            end
                        
            if isAmp
                if Nta ~= length(beta)
                    err = 1;
                end
                if err
                    error('waves must be of size Nsig x Nwg x Nbeta');
                end
            else
                if length(time) ~= Nta
                    err = 1;
                end
                if err
                    error('waves must be of size Nsig x Nwg x Ntime');
                end
            end
            
            tda.wgPos = wgPos;
            
            if isAmp
                tda.waveBeta = beta;
                tda.waveAmps = cell(Nsig, Nwg);
                for m = 1:Nsig
                    for n = 1:Nwg
                        tda.waveAmps{m, n} = squeeze(waves(m, n, :));
                    end
                end
            else               
                tda.waveTime = time;
                tda.waveSigs = cell(Nsig, Nwg);
                for m = 1:Nsig
                    for n = 1:Nwg
                        tda.waveSigs{m, n} = squeeze(waves(m, n, :));
                    end
                end
            end
        end
        
        function [] = SetMotionTimeLimits(tda, sigInds, dofs, startTime, stopTime)
            if isempty(tda.motTimeLims)
                tda.motTimeLims = cell(tda.nSig, tda.motDof);
            end
            for m = 1:length(sigInds)
                for n = 1:length(dofs)
                    tda.motTimeLims{sigInds(m), dofs(n)} = [startTime stopTime];
                end
            end
            tda.smotions = [];
        end
        
        function [] = SetPtoTimeLimits(tda, sigInds, dofs, startTime, stopTime)
            if isempty(tda.ptoTimeLims)
                tda.ptoTimeLims = cell(tda.nSig, tda.ptoDof);
            end
            for m = 1:length(sigInds)
                for n = 1:length(dofs)
                    tda.ptoTimeLims{sigInds(m), dofs(n)} = [startTime stopTime];
                end
            end
            tda.sptoKinematic = [];
            tda.sptoDynamic = [];
        end
        
        function [] = SetWaveTimeLimits(tda, startTime, stopTime)
            tda.waveTimeLims = [startTime stopTime];
            tda.swaveSigs = [];
        end
        
        function [mots, timeFreq] = GetMotions(tda, sigInds, dofs, varargin)
            [mots, timeFreq] = tda.getValues('motions', sigInds, dofs, varargin{:});
        end
        
        function [ptoKin, timeFreq] = GetPtoKinematic(tda, sigInds, dofs, varargin)
            [ptoKin, timeFreq] = tda.getValues('ptoKinematic', sigInds, dofs, varargin{:});
        end
        
        function [ptoDyn, timeFreq] = GetPtoDynamic(tda, sigInds, dofs, varargin)
            [ptoDyn, timeFreq] = tda.getValues('ptoDynamic', sigInds, dofs, varargin{:});
        end
        
        function [pow, time] = Power(tda, sigInds, dofs, varargin)
            [pow, time] = tda.power(sigInds, dofs, varargin{:});
        end
        
        function [waves, timeFreq, wgPos, beta] = GetWaves(tda, sigInds, wgInds, varargin)
            [opts, args] = checkOptions({{'filter'}, {'amps'}, {'spectra'}, {'fullTime'}}, varargin);
            
            if isempty(tda.wgPos)
                waves = [];
                timeFreq = [];
                wgPos = [];
                beta = [];
                return;
            end
            
            filt = opts(1);
            isAmps = opts(2);
            spec = opts(3);
            fullTime = opts(4);
            
            if ischar(sigInds)
                if strcmp(sigInds, ':');
                    sigInds = 1:tda.nSig;
                end
            end
            
            if ischar(wgInds)
                if strcmp(wgInds, ':');
                    [Nwg, ~] = size(tda.wgPos);
                    wgInds = 1:Nwg;
                end
            end
            
            Nsig = length(sigInds);
            Nwg = length(wgInds);
                        
            waves = cell(Nsig, Nwg);
            beta = tda.waveBeta;
            
            if filt && isempty(tda.fwaveSigs)
                tda.computeFilter;
            end
            
            if isAmps && isempty(tda.waveAmps)
                tda.computeWaveAmps;
            end
            
            if spec && isempty(tda.specWaves)
                tda.computeSpec;
            end
            
            tlims = tda.waveTimeLims;
            
            if spec
                timeFreq = tda.waveFreq;
            else
                timeFreq = tda.waveTime;
                [iStart, iStop] = tda.tlimInds(tlims, timeFreq, fullTime);
                timeFreq = timeFreq(iStart:iStop);
            end
                        
            wgPos = zeros(length(wgInds), 2);
            for m = 1:Nsig
                for n = 1:Nwg
                    if isempty(tda.waveSigs) || isAmps
                        waves{m, n} = tda.waveAmps{sigInds(m), wgInds(n)};
                    else
                        if filt
                            sigmn = tda.fwaveSigs{sigInds(m), wgInds(n)};
                            waves{m, n} = sigmn(iStart:iStop);
                        elseif spec
                            waves{m, n} = tda.swaveSigs{sigInds(m), wgInds(n)};
                        else
                            sigmn = tda.waveSigs{sigInds(m), wgInds(n)};
                            waves{m, n} = sigmn(iStart:iStop);
                        end
                    end
                    wgPos(n, :) = tda.wgPos(wgInds(n), :);
                end
            end
        end
        
        function [startTime, stopTime] = GetMotionTimeLimits(tda, sigInd, dof)
            startTime = tda.motTimeLims{sigInd, dof}(1);
            stopTime = tda.motTimeLims{sigInd, dof}(2);
        end
        
        function [startTime, stopTime] = GetPtoTimeLimits(tda, sigInd, dof)
            startTime = tda.ptoTimeLims{sigInd, dof}(1);
            stopTime = tda.ptoTimeLims{sigInd, dof}(2);
        end
        
        function [startTime, stopTime] = GetWaveTimeLimits(tda)
            startTime = tda.waveTimeLims(1);
            stopTime = tda.waveTimeLims(2);
        end
    end
    
    methods (Access = protected)
        function [] = setValues(tda, type, dofs, time, signals)
            
            switch type
                case 'motions'
                    sType = 'mot';
                    sigName = 'Motions';
                case 'ptoKinematic'
                    sType = 'pto';
                    sigName = 'PtoKinematic';
                case 'ptoDynamic'
                    sType = 'pto';
                    sigName = 'PtoDynamic';
                otherwise
                    error('setValues type not recognized');
            end
            
            if length(dofs) == 1 && ndims(signals) == 2
                [Nsig, Nt] = size(signals);
                sigsTemp = signals;
                signals = zeros(Nsig, 1, Nt);
                signals(:, 1, :) = sigsTemp;
            elseif length(dofs) > 1 && ndims(signals) == 2
               [Ndof, Nt] = size(signals);
               sigsTemp = signals;
               signals = zeros(1, Ndof, Nt);
               signals(1, :, :) = sigsTemp; 
            end
            
            [Nsig, Ndof, Nt] = size(signals);
            err = 0;
            if isempty(tda.nSig)
                tda.nSig = Nsig;
            elseif (tda.nSig ~= Nsig)
                err = 1;
            end
            
            if length(dofs) ~= Ndof
                err = 1;
            end
            
            if eval(['isempty(tda.' sType 'Dof) || tda.' sType 'Dof < max(dofs)'])
                eval(['tda.' sType 'Dof = max(dofs);']);
            end
            
            if length(time) ~= Nt
                err = 1;
            end
            
            if err
                error([sigName ' must be of size Nsig x Ndof x Ntime']);
            end
            
            if eval(['isempty(tda.' sType 'Time)'])
                eval(['tda.' sType 'Time = cell(tda.nSig, tda.' sType 'Dof);'])
            end
            
            if eval(['isempty(tda.' type ')'])
                eval(['tda.' type ' = cell(tda.nSig, tda.' sType 'Dof);'])
            end
            
            for m = 1:tda.nSig
                for n = 1:Ndof
                    eval(['tda.' sType 'Time{m, dofs(n)} = time;'])
                    eval(['tda.' type '{m, dofs(n)} = squeeze(signals(m,n,:));'])
                end
            end
        end
        
        function [sigs, timeFreq] = getValues(tda, type, sigInds, dofs, varargin)
            
            [opts, args] = checkOptions({{'filter'}, {'spectra'}, {'fullTime'}}, varargin);
            
            filt = opts(1);
            spec = opts(2);
            fullTime = opts(3);
                        
            fstr = '';
            if filt
                fstr = 'f';
            elseif spec
                fstr = 's';
            end
            
            switch type
                case 'motions'
                    sType = 'mot';
                case 'ptoKinematic'
                    sType = 'pto';
                case 'ptoDynamic'
                    sType = 'pto';
                otherwise
                    error('getValues type not recognized');
            end
            
            if ischar(sigInds)
                if strcmp(sigInds, ':');
                    sigInds = 1:tda.nSig;
                end
            end
            
            if ischar(dofs)
                if strcmp(dofs, ':');
                    eval(['dofs = 1:tda.' sType 'Dof;'])
                end
            end
            
            if max(sigInds) > tda.nSig
                error('The signal indices exceed the number of signals');
            end
            
            if eval(['max(dofs) > tda.' sType 'Dof'])
                error(['The dof indices exceed the number of ' sType ' dof']);
            end
            
            if eval(['filt && isempty(tda.' fstr type ')'])
                tda.computeFilter;
            end
            
            if eval(['spec && isempty(tda.' fstr type ')'])
                tda.computeSpec;
            end
            
            Nsig = length(sigInds);
            Ndof = length(dofs);
            sigs = cell(Nsig, Ndof);
            timeFreq = cell(Nsig, Ndof);
            
            eval(['tlims = tda.' sType 'TimeLims;']);
            
            if isempty(tlims)
                tlims = cell(Nsig, Ndof);
            end
            
            for m = 1:Nsig
                for n = 1:Ndof
                    if spec
                        eval(['timeFreq{m, n} = tda.' sType 'Freq{sigInds(m), dofs(n)};'])
                        eval(['sigs{m, n} = tda.' fstr type '{sigInds(m), dofs(n)};'])
                    else
                        eval(['timemn = tda.' sType 'Time{sigInds(m), dofs(n)};'])
                        eval(['sigmn = tda.' fstr type '{sigInds(m), dofs(n)};'])
                        
                        [iStart, iStop] = tda.tlimInds(tlims{m, n}, timemn, fullTime);
                        
                        timeFreq{m, n} = timemn(iStart:iStop);
                        sigs{m, n} = sigmn(iStart:iStop);
                    end  
                end
            end
        end
        
        function [pow, time] = power(tda, sigInds, dofs, varargin)
            [ptoKin] = tda.getValues('ptoKinematic', sigInds, dofs, varargin{:});
            [ptoDyn, time] = tda.getValues('ptoDynamic', sigInds, dofs, varargin{:});
            
            [Nsig, Ndof] = size(ptoKin);
            pow = cell(Nsig, Ndof);
            
            for m = 1:Nsig
                for n = 1:Ndof
                    pow{m, n} = ptoKin{m, n}.*ptoDyn{m, n};
                end
            end
        end
        
        function [] = computeFilter(tda)
        end
        
        function [] = computeSpec(tda)
            
            if ~isempty(tda.motions)
                tda.smotions = cell(tda.nSig, tda.motDof);
                tda.motFreq = cell(tda.nSig, tda.motDof);
                tlims = tda.motTimeLims;
                if isempty(tlims)
                    tlims = cell(tda.nSig, tda.motDof);
                end
                for m = 1:tda.nSig
                    for n = 1:tda.motDof
                        timemn = tda.motTime{m, n};
                        sigmn = tda.motions{m, n};
                        [iStart, iStop] = tda.tlimInds(tlims{m, n}, timemn, false);
                        
                        [spec, freq] = TimeDomainAnalysis.FFT(timemn(iStart:iStop), sigmn(iStart:iStop));
                        tda.smotions{m, n} = spec;
                        tda.motFreq{m, n} = freq;
                    end
                end
            end
            
            if ~isempty(tda.ptoKinematic)
                tda.sptoKinematic = cell(tda.nSig, tda.ptoDof);
                tda.ptoFreq = cell(tda.nSig, tda.ptoDof);
                tlims = tda.ptoTimeLims;
                if isempty(tlims)
                    tlims = cell(tda.nSig, tda.ptoDof);
                end
                for m = 1:tda.nSig
                    for n = 1:tda.ptoDof
                        timemn = tda.ptoTime{m, n};
                        sigmn = tda.ptoKinematic{m, n};
                        [iStart, iStop] = tda.tlimInds(tlims{m,n}, timemn, false);
                        
                        [spec, freq] = TimeDomainAnalysis.FFT(timemn(iStart:iStop), sigmn(iStart:iStop));
                        tda.sptoKinematic{m, n} = spec;
                        tda.ptoFreq{m, n} = freq;
                    end
                end
            end
            
            if ~isempty(tda.ptoDynamic)
                tda.sptoDynamic = cell(tda.nSig, tda.ptoDof);
                tda.ptoFreq = cell(tda.nSig, tda.ptoDof);
                tlims = tda.ptoTimeLims;
                if isempty(tlims)
                    tlims = cell(tda.nSig, tda.ptoDof);
                end
                for m = 1:tda.nSig
                    for n = 1:tda.ptoDof
                        timemn = tda.ptoTime{m, n};
                        sigmn = tda.ptoDynamic{m, n};
                        [iStart, iStop] = tda.tlimInds(tlims{m,n}, timemn, false);
                        
                        [spec, freq] = TimeDomainAnalysis.FFT(timemn(iStart:iStop), sigmn(iStart:iStop));
                        tda.sptoDynamic{m, n} = spec;
                        tda.ptoFreq{m, n} = freq;
                    end
                end
            end
            
            if ~isempty(tda.waveSigs)
                [Nwg, ~] = size(tda.wgPos);
                tda.swaveSigs = cell(tda.nSig, Nwg);
                tda.waveFreq = cell(tda.nSig, Nwg);
                for m = 1:tda.nSig
                    for n = 1:Nwg
                        timemn = tda.waveTime{m, n};
                        sigmn = tda.waveSigs{m, n};
                        [iStart, iStop] = tda.tlimInds(tda.waveTimeLims, timemn, false);
                        
                        [spec, freq] = TimeDomainAnalysis.FFT(timemn(iStart:iStop), sigmn(iStart:iStop));
                        tda.swaveSigs{m, n} = spec;
                        tda.waveFreq{m, n} = freq;
                    end
                end
            end
        end
        
        function [] = computeWaveAmps(tda)
        end
        
        function [iStart, iStop] = tlimInds(tda, tlims, time, fullTime)
            if isempty(tlims) || fullTime
                iStart = 1;
                iStop = length(time);
            else
                if isinf(tlims(1))
                    iStart = 1;
                else
                    iStart = indexOf(time, tlims(1));
                end
                if isinf(tlims(2))
                    iStop = length(time);
                else
                    iStop = indexOf(time, tlims(2));
                end
            end
        end
    end
    
    methods (Static)
        function [spec, freq, offset] = FFT(time, signal)
            N = length(signal);

            offset = mean(signal);

            S = fft(signal)/N;
            spec = 2*S(1:floor(N/2)+1);
            dt = time(2)-time(1);
            t0 = time(1);
            freq = ((0:floor(N/2))/dt/N).';
            % phase shift base on first time value
            spec = spec.*exp(1i*2*pi*freq*t0);
        end
    end
end