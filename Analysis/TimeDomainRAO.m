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
classdef TimeDomainRAO < TimeDomainAnalysis
    
    properties (Access = private)
        waveAmps;
        waveBeta;
        freqs;
    end
    
    properties (Dependent)
        Frequencies;
    end
    
    methods
        
        function [tda] = TimeDomainRAO()
        end
        
        function [] = set.Frequencies(tda, val)
            Nf = length(val);
            if isempty(tda.nSig)
                tda.nSig = Nf;
            else
                if tda.nSig ~= Nf
                    error('The number of frequencies must be equal to the number of signals.');
                end
            end
            tda.freqs = val;
        end
        function [val] = get.Frequencies(tda)
            val = tda.freqs;
        end
                
        function [mots, timeFreq] = GetMotions(tda, sigInds, dofs, varargin)
            [mots, timeFreq] = tda.raoGetValues('motions', sigInds, dofs, varargin{:});
        end

        function [hinLoads, timeFreq] = GetHingeLoads(tda, sigInds, dofs, varargin)
            [hinLoads, timeFreq] = tda.raoGetValues('hingeLoads', sigInds, dofs, varargin{:});
        end
        
        function [ptoKin, timeFreq] = GetPtoKinematic(tda, sigInds, dofs, varargin)
            [ptoKin, timeFreq] = tda.raoGetValues('ptoKinematic', sigInds, dofs, varargin{:});
        end
        
        function [ptoDyn, timeFreq] = GetPtoDynamic(tda, sigInds, dofs, varargin)
            [ptoDyn, timeFreq] = tda.raoGetValues('ptoDynamic', sigInds, dofs, varargin{:});
        end

        function [ptoPow, timeFreq] = GetPtoPower(tda, sigInds, dofs, varargin)
            [ptoPow, timeFreq] = tda.raoGetValues('ptoPower', sigInds, dofs, varargin{:});
        end
        
        function [pow, timeFreq] = Power(tda, sigInds, dofs, varargin)
            [opts, args] = checkOptions({{'rao'}, {'noNorm'}, {'CW', 1}, {'mat'}, {'mean'}}, varargin);
            
            rao = opts(1);
            noNorm = opts(2);
            compCW = opts(3);
            mat = opts(4);
            if opts(5)
                rao = true;
                noNorm = true;
            end
            
            ai = tda.GetWaves(sigInds, 1, 'amps', 'mat');
            if isempty(ai)
                noNorm = true;
            end
            
            if ~noNorm
                [~, Nbeta] = size(ai);
                if Nbeta > 1
                    ai = ai(:,1);
                end
                norm = abs(ai).^2;
                if compCW
                    rho = args{3};
                    waves = PlaneWaves(abs(ai), 1./tda.freqs, 0, tda.h);
                    norm = waves.EnergyFlux(rho);
                end
            else
                norm = ones(size(tda.freqs));
            end
            
            [pow0, timeFreq] = tda.power(sigInds, dofs, varargin{:});
            
            if rao
                timeFreq = tda.freqs;
                [Nf, Ndof] = size(pow0);
                pow = cell(Ndof, 1);
                for m = 1:Ndof
                    pow{m} = zeros(1, Nf);
                    for n = 1:Nf
                        pow{m}(n) = mean(pow0{n, m})/norm(n);
                    end
                end
            else
                pow = pow0;
            end
            
            if mat
                if rao
                    powM = zeros(Ndof, Nf);
                    for m = 1:Ndof
                        powM(m,:) = pow{m};
                    end
                else
                    Nt = length(pow{1,1});
                    powM = zeros(Nf, Ndof, Nt);
                    for m = 1:Nf
                        for n = 1:Ndof
                            powM(m,n,:) = pow{m,n};
                        end
                    end
                end
                pow = powM;
            end
        end
        
        function [waves, timeFreq, wgPos, beta] = GetWaves(tda, sigInds, wgInds, varargin)
            [opts, args] = checkOptions({{'amps'}, {'flux', 1}, {'betaInd', 1}, {'mat'}}, varargin);
            
            if isempty(tda.wgPos)
                waves = [];
                timeFreq = [];
                wgPos = [];
                beta = [];
                return;
            end
            
            isAmps = opts(1);
            rho = [];
            T = 1./tda.freqs;
            if opts(2)
                rho = args{2};
            end
            betaInds = 1:length(tda.waveBeta);
            if opts(3)
                betaInds = args{3};
            end
            mat = opts(4);
            
            if isempty(tda.waveSigs)
                isAmps = true;
            end
            
            if isAmps
                [sigInds, wgInds, Nsig, Nwg] = tda.getWaveInds(sigInds, wgInds);
            
                if isempty(tda.waveAmps)
                    tda.computeWaveAmps;
                end
            
                timeFreq = [];
                beta = tda.waveBeta;
                waves = cell(Nsig, Nwg);
                wgPos = zeros(length(wgInds), 2);
                
                k0 = IWaves.SolveForK(2*pi*tda.freqs, tda.h);
                for m = 1:Nsig
                    for n = 1:Nwg
                        
                        wgPos(n, :) = tda.wgPos(wgInds(n), :);
                        r0 = tda.meanPos - wgPos(n, :);
                        amps0 = tda.waveAmps{sigInds(m), wgInds(n)};
                        
                        % adjust phase to WEC position
                        amps = zeros(1,length(betaInds));
                        for o = 1:length(betaInds)
                            k = k0(n)*[cos(beta(betaInds(o))), sin(beta(betaInds(o)))];
                            amps(o) = amps0(o)*exp(-1i*dot(k, r0));
                        end

                        % energy flux
                        if ~isempty(rho)
                            wav = PlaneWaves(abs(amps(betaInds)), T(m)*ones(size(amps(betaInds))), tda.waveBeta(betaInds), tda.h);
                            amps = wav.EnergyFlux(rho)';
                        end
                        waves{m, n} = amps;
                    end
                end
                
                if mat
                    Nb = length(waves{1,1});
                    wavesM = zeros(Nsig, Nwg, Nb);
                    for m = 1:Nsig
                        for n = 1:Nwg
                            wavesM(m,n,:) = waves{m,n};
                        end
                    end
                    waves = wavesM;
                    waves = squeeze(waves);
                end
            else
                beta = [];
                [waves, timeFreq, wgPos] = tda.getWaves(sigInds, wgInds, varargin{:});
            end
        end
        
        function [] = SetWaves(tda, wgPos, time, waves, varargin)
            isAmp = false;
            if isempty(time)
                isAmp = true;
            end
            
            if isAmp
                if isempty(varargin)
                    beta = 0;
                else
                    beta = varargin{1};
                end
            end
            
            [waves, Nsig, Nwg, Nta] = tda.checkWaves(wgPos, time, waves);
            
            if isAmp
                if Nta ~= length(beta)
                    error('waves must be of size Nsig x Nwg x Nbeta');
                end
                
                tda.wgPos = wgPos;
                tda.waveBeta = beta;
                tda.waveAmps = cell(Nsig, Nwg);
                                
                for m = 1:Nsig
                    for n = 1:Nwg
                        tda.waveAmps{m, n} = squeeze(waves(m, n, :));
                    end
                end
            else
                tda.setTimeWaves(wgPos, time, waves, Nsig, Nwg, waveRamp);
            end
        end
    
    end
    
    methods (Access = protected)
        function [sigs, timeFreq] = raoGetValues(tda, type, sigInds, dofs, varargin)
            
            [opts, args] = checkOptions({{'rao', 1}, {'noNorm'}, {'normPhase'}, ...
                {'offset', 2}, {'mat'}, {'normPhaseSurge'}}, varargin);
            
            rao = opts(1);
            noNorm = opts(2);
            normPhase = opts(3);
            offset = false;
            if opts(4)
                offset = true;
                beginTime = [];
                if args{4}{1}
                    beginTime = args{4}{2};
                end
            end
            mat = opts(5);
            normSur = opts(6);
            
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
                    switch type
                        case 'hingeLoads'
                            dofs = 1:tda.forceDof;
                        case 'motions'
                            dofs = 1:tda.motDof;
                    end
                end
            end
            
            ai = tda.GetWaves(sigInds, 1, 'amps', 'betaInd', 1, 'mat');
            
            if isempty(ai)
                noNorm = true;
            end
            
            if noNorm
                ai = ones(size(tda.freqs));
            elseif normPhase
                ai = exp(1i*angle(ai));
            end
            
            if rao
                [specs, freq] = tda.getValues(type, sigInds, dofs, 'spectra');
                Na = args{1};
                if ~isfloat(Na)
                    Na = 1;
                elseif Na < 1
                    Na = 1;
                end
                [Nf, Ndof] = size(specs);
                sigs = cell(Ndof, Na);
                timeFreq = tda.freqs;
                
                if normSur
                    beta = tda.waveBeta(1);
                    k0 = IWaves.SolveForK(2*pi*tda.freqs, tda.h);
                    for m = 1:Nf
                        k = k0(m)*[cos(beta), sin(beta)];
                        xpos = specs{m,1}(1);
                        r0 = [xpos, 0];
                        ai(m) = ai(m)*exp(-1i*dot(k, r0));
                    end
                end
                
                if offset
                    Na = 1;
                    [fullTime, time] = tda.getValues(type, sigInds, dofs, 'fullTime');
                end
                                
                for m = 1:Ndof
                    for n = 1:Na
                        sigs{m, n} = zeros(1, Nf);
                        for o = 1:Nf
                            if offset
                                delta = 0;
                                if ~isempty(beginTime)
                                    iTime = [indexOf(time{o, m}, beginTime(1)), indexOf(time{o, m}, beginTime(2))];
                                    if iTime(2) ~= iTime(1)
                                        delta = mean(fullTime{o, m}(iTime(1):iTime(2)));
                                    end
                                end
                                sigs{m, 1}(o) = (specs{o, m}(1) - delta)/ai(o);
                                sigs{m, 2}(o) = delta/ai(o);
                            else
                                indf = indexOf(freq{o,m}, n*tda.freqs(o));
                                sigs{m, n}(o) = specs{o, m}(indf)/ai(o);
                            end
                        end
                    end
                end
                
                if mat
                    sigsM = zeros(Ndof, Na, Nf);
                    for m = 1:Ndof
                        for n = 1:Na
                            for o = 1:Nf
                                sigsM(m, n, o) = sigs{m, n}(o);
                            end
                        end
                    end
                    sigs = squeeze(sigsM);
                end
            else
                [sigs, timeFreq] = tda.getValues(type, sigInds, dofs, varargin{:});
            end
            
        end
        
        function [] = computeWaveAmps(tda)
        end
    end
end