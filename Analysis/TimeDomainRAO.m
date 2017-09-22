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
        freqs;
        h;
    end
    
    properties (Dependent)
        Frequencies;
        WaterDepth;
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
            val = ras.freqs;
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
        
        function [mots, timeFreq] = GetMotions(tda, sigInds, dofs, varargin)
            [mots, timeFreq] = tda.raoGetValues('motions', sigInds, dofs, varargin{:});
        end
        
        function [ptoKin, timeFreq] = GetPtoKinematic(tda, sigInds, dofs, varargin)
            [ptoKin, timeFreq] = tda.raoGetValues('ptoKinematic', sigInds, dofs, varargin{:});
        end
        
        function [ptoDyn, timeFreq] = GetPtoDynamic(tda, sigInds, dofs, varargin)
            [ptoDyn, timeFreq] = tda.raoGetValues('ptoDynamic', sigInds, dofs, varargin{:});
        end
        
        function [pow, timeFreq] = Power(tda, sigInds, dofs, varargin)
            [opts, args] = checkOptions({{'rao'}, {'noNorm'}, {'CW', 1}}, varargin);
            
            noNorm = opts(2);
            compCW = opts(3);
            
            ai = tda.GetWaves(sigInds, 1, 'amps');
            if isempty(ai)
                noNorm = true;
            end
            
            if ~noNorm
                ai = [ai{:}].';
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
            
            if opts(1)
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
        end
    
    end
    
    methods (Access = protected)
        function [sigs, timeFreq] = raoGetValues(tda, type, sigInds, dofs, varargin)
            
            [opts, args] = checkOptions({{'rao', 1}, {'noNorm'}}, varargin);
            
            noNorm = opts(2);
            [ai, ~, wgPos, beta] = tda.GetWaves(sigInds, 1, 'amps');
            
            if isempty(ai)
                noNorm = true;
            end
            
            if ~noNorm
                ai = [ai{:}].';
                r0 = tda.meanPos - wgPos;
                k0 = IWaves.SolveForK(2*pi*tda.freqs, tda.h);
                k = k0*[cos(beta), sin(beta)];
                Nf = length(tda.freqs);
                for n = 1:Nf
                    ai(n) = ai(n)*exp(-1i*dot(k(n,:), r0));
                end
            else
                ai = ones(size(tda.freqs));
            end
            
            if opts(1)
                [specs, freq] = tda.getValues(type, sigInds, dofs, 'spectra');
                Na = args{1};
                [Nf, Ndof] = size(specs);
                sigs = cell(Ndof, Na);
                timeFreq = tda.freqs;
                                
                for m = 1:Ndof
                    for n = 1:Na
                        sigs{m, n} = zeros(1, Nf);
                        for o = 1:Nf
                            indf = indexOf(freq{o,m}, n*tda.freqs(o));
                            sigs{m, n}(o) = specs{o, m}(indf)/ai(o);
                        end
                    end
                end
            else
                [sigs, timeFreq] = tda.getValues(type, sigInds, dofs, varargin{:});
            end
        end
    end
end