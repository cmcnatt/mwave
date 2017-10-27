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
        sigDes;
        motDof;
        ptoDof;
        motions;
        ptoKinematic;
        ptoDynamic;
        waveSigs;
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
        meanPos;
        h;
    end
    
    properties (Dependent)
        MotionDoF;
        PtoDoF;
        SignalDescription;
        NSignals;
        MeanBodyPosition;
        WaterDepth;
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
        
        function [] = set.SignalDescription(tda, val)
            tda.sigDes = val;
        end
        function [val] = get.SignalDescription(tda)
            val = tda.sigDes;
        end
        
        function [] = set.MeanBodyPosition(tda, val)
            tda.meanPos = val;
        end
        function [val] = get.MeanBodyPosition(tda)
            val = tda.meanPos;
        end
        
        function [] = set.WaterDepth(tda, val)
            if ~isscalar(val)
                error('WaterDepth must be a scalar');            
            end
            tda.h = val;
        end
        function [val] = get.WaterDepth(tda)
            val = tda.h;
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
            
            [waves, Nsig, Nwg] = tda.checkWaves(wgPos, time, waves, varargin{:});
            
            tda.setTimeWaves(wgPos, time, waves, Nsig, Nwg);
        end
        
        function [] = SetMotionTimeLimits(tda, sigInds, dofs, startTime, stopTime)
            if isempty(tda.motTimeLims)
                tda.motTimeLims = cell(tda.nSig, tda.motDof);
            end
            
            if isempty(sigInds)
                sigInds = 1:tda.nSig;
            elseif ischar(sigInds)
                if strcmp(sigInds, ':');
                    sigInds = 1:tda.nSig;
                end
            end
            
            if isempty(dofs)
                dofs = 1:tda.motDof;
            elseif ischar(dofs)
                if strcmp(dofs, ':');
                    dofs = 1:tda.motDof;
                end
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
            
            if isempty(sigInds)
                sigInds = 1:tda.nSig;
            elseif ischar(sigInds)
                if strcmp(sigInds, ':');
                    sigInds = 1:tda.nSig;
                end
            end
            
            if isempty(dofs)
                dofs = 1:tda.ptoDof;
            elseif ischar(dofs)
                if strcmp(dofs, ':');
                    dofs = 1:tda.ptoDof;
                end
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
            [opts] = checkOptions({{'mean'}}, varargin);
            
            [pow, time] = tda.power(sigInds, dofs, varargin{:});
            
            if opts(1)
                [M, N] = size(pow);
                meanPow = zeros(M, N);
                for m = 1:M
                    for n = 1:N
                        meanPow(m, n) = mean(pow{m, n});
                    end
                end
                pow = meanPow;
            end
        end
        
        function [damping] = GetPtoDamping(tda, sigInds, dofs, varargin)
            [ptoDyn] = tda.getValues('ptoDynamic', sigInds, dofs, varargin{:});
            [ptoKin] = tda.getValues('ptoKinematic', sigInds, dofs, varargin{:});
            
            [Nsig, Ndof] = size(ptoDyn);
            damping = zeros(Nsig, Ndof);
            
            for m = 1:Nsig
                for n = 1:Ndof
                    damping(m, n) = sum(ptoDyn{m,n}.*ptoKin{m,n})./sum(ptoKin{m,n}.^2);
                end
            end
        end
        
        function [waves, timeFreq, wgPos] = GetWaves(tda, sigInds, wgInds, varargin)
            [waves, timeFreq, wgPos] = tda.getWaves(sigInds, wgInds, varargin{:});
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
            
            [opts, args] = checkOptions({{'filter'}, {'spectra'}, {'fullTime'}, {'noMean'}}, varargin);
            
            filt = opts(1);
            spec = opts(2);
            fullTime = opts(3);
            noMean = opts(4);
                        
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
            
            if isempty(sigInds)
                sigInds = 1:tda.nSig;
            elseif ischar(sigInds)
                if strcmp(sigInds, ':');
                    sigInds = 1:tda.nSig;
                end
            end
            
            if isempty(dofs)
                eval(['dofs = 1:tda.' sType 'Dof;'])
            elseif ischar(dofs)
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
            
            if spec
                tda.computeSpec(varargin{:});
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
        
        function [pow, timeFreq] = power(tda, sigInds, dofs, varargin)
            [opts, args] = checkOptions({{'spectra'}, {'smooth', 1}}, varargin);
            
            spectra = opts(1);
            smooth = opts(2);
            windowDf = [];
            if smooth
                windowDf = args{2};
                for n = 1:length(varargin)
                    if strcmp(varargin{n}, 'smooth')
                        varargin{n} = []; 
                        varargin{n+1} = []; 
                        break; 
                    end
                end
            end

            [ptoKin] = tda.getValues('ptoKinematic', sigInds, dofs, varargin{:});
            [ptoDyn, timeFreq] = tda.getValues('ptoDynamic', sigInds, dofs, varargin{:});
            
            [Nsig, Ndof] = size(ptoKin);
            pow = cell(Nsig, Ndof);
            
            for m = 1:Nsig
                for n = 1:Ndof
                    if spectra
                        spec = real(0.5*ptoKin{m, n}.*conj(ptoDyn{m, n}));
                        if smooth
                            spec = TimeDomainAnalysis.SmoothSpectrum(timeFreq{m,n}, spec, windowDf);
                        end
                        pow{m, n} = real(spec);
                    else
                        pow{m, n} = ptoKin{m, n}.*ptoDyn{m, n};
                    end
                end
            end
        end
        
        function [waves, timeFreq, wgPos] = getWaves(tda, sigInds, wgInds, varargin)
            [opts, args] = checkOptions({{'filter'}, {'spectra'}, {'fullTime'}, {'mat'}}, varargin);
            
            if isempty(tda.wgPos)
                waves = [];
                timeFreq = [];
                wgPos = [];
                return;
            end
            
            filt = opts(1);
            spec = opts(2);
            fullTime = opts(3);
            mat = opts(4);
            
            [sigInds, wgInds, Nsig, Nwg] = tda.getWaveInds(sigInds, wgInds);
                        
            waves = cell(Nsig, Nwg);
            
            if filt && isempty(tda.fwaveSigs)
                tda.computeFilter;
            end
                        
            if spec && isempty(tda.specWaves)
                tda.computeSpec(varargin{:});
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
                    if filt
                        sigmn = tda.fwaveSigs{sigInds(m), wgInds(n)};
                        waves{m, n} = sigmn(iStart:iStop);
                    elseif spec
                        waves{m, n} = tda.swaveSigs{sigInds(m), wgInds(n)};
                    else
                        sigmn = tda.waveSigs{sigInds(m), wgInds(n)};
                        waves{m, n} = sigmn(iStart:iStop);
                    end
                    wgPos(n, :) = tda.wgPos(wgInds(n), :);
                end
            end
            
            if mat
                Nt = length(waves{1,1});
                wavesM = zeros(Nsig, Nwg, Nt);
                for m = 1:Nsig
                    for n = 1:Nwg
                        wavesM(m,n,:) = waves{m,n};
                    end
                end
                waves = wavesM;
            end
        end
        
        function [sigInds, wgInds, Nsig, Nwg] = getWaveInds(tda, sigInds, wgInds)
            if isempty(sigInds)
                sigInds = 1:tda.nSig;
            elseif ischar(sigInds)
                if strcmp(sigInds, ':');
                    sigInds = 1:tda.nSig;
                end
            end
            
            [Nwg, ~] = size(tda.wgPos);
            if isempty(wgInds)
                wgInds = 1:Nwg;
            elseif ischar(wgInds)
                if strcmp(wgInds, ':');
                    wgInds = 1:Nwg;
                end
            end
            
            Nsig = length(sigInds);
            Nwg = length(wgInds);
        end
        
        function [waves, Nsig, Nwg, Nta] = checkWaves(tda, wgPos, time, waves, varargin)
            [Nwg, ~] = size(wgPos);
                        
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
                        
            if ~isempty(time)
                if length(time) ~= Nta
                    err = 1;
                end
                if err
                    error('waves must be of size Nsig x Nwg x Ntime');
                end
            end
        end
        
        function [] = setTimeWaves(tda, wgPos, time, waves, Nsig, Nwg)
            
            tda.wgPos = wgPos;
            
            tda.waveTime = time;
            tda.waveSigs = cell(Nsig, Nwg);
            for m = 1:Nsig
                for n = 1:Nwg
                    tda.waveSigs{m, n} = squeeze(waves(m, n, :));
                end
            end
        end
        
        function [] = computeFilter(tda)
        end
        
        function [] = computeSpec(tda, varargin)
            
            [opts, args] = checkOptions({{'noMean'}, {'smooth', 1}}, varargin);
            
            noMean = opts(1);
            windowDf = [];
            if opts(2)
                windowDf = args{2};
            end
            
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
                        
                        [spec, freq] = TimeDomainAnalysis.FFT(timemn(iStart:iStop), sigmn(iStart:iStop), noMean);
                        if ~isempty(windowDf)
                            spec = TimeDomainAnalysis.SmoothSpectrum(freq, spec, windowDf);
                        end
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
                        
                        [spec, freq] = TimeDomainAnalysis.FFT(timemn(iStart:iStop), sigmn(iStart:iStop), noMean);
                        if ~isempty(windowDf)
                            spec = TimeDomainAnalysis.SmoothSpectrum(freq, spec, windowDf);
                        end
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
                        
                        [spec, freq] = TimeDomainAnalysis.FFT(timemn(iStart:iStop), sigmn(iStart:iStop), noMean);
                        if ~isempty(windowDf)
                            spec = TimeDomainAnalysis.SmoothSpectrum(freq, spec, windowDf);
                        end
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
                        
                        [spec, freq] = TimeDomainAnalysis.FFT(timemn(iStart:iStop), sigmn(iStart:iStop), noMean);
                        if ~isempty(windowDf)
                            spec = TimeDomainAnalysis.SmoothSpectrum(freq, spec, windowDf);
                        end
                        tda.swaveSigs{m, n} = spec;
                        tda.waveFreq{m, n} = freq;
                    end
                end
            end
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
        function [spec, freq, offset] = FFT(time, signal, varargin)
            noMean = false;
            if ~isempty(varargin)
                noMean = varargin{1};
            end
            
            N = length(signal);

            offset = mean(signal);
            if noMean
                signal = signal - offset;
            end

            S = fft(signal)/N;
            spec = 2*S(1:floor(N/2)+1);
            spec(1) = spec(1)/2;
            dt = time(2)-time(1);
            t0 = time(1);
            freq = ((0:floor(N/2))/dt/N).';
            % phase shift base on first time value
            spec = spec.*exp(1i*2*pi*freq*t0);
        end
        
        function [specOut] = SmoothSpectrum(freq, spectrum, windowDf)
            winN = indexOf(freq-freq(1), windowDf);
            
            if mod(winN,2) == 0
                winN = winN + 1;
            end
            delN = (winN-1)/2;
            
            specOut = zeros(size(spectrum));
            N = length(spectrum);
            for n = 1:N
                iStart = n - delN;
                iStop = n + delN;
                if iStart < 1
                    iStart = 1;
                    iStop = n + (n - iStart);
                end
                if iStop > N
                    iStop = N;
                    iStart = n - (n - iStart) + 1;
                end
                specSeg = spectrum(iStart:iStop);
                specOut(n) = mean(abs(specSeg))*exp(1i*angle(spectrum(n)));
            end
                
        end
    end
end