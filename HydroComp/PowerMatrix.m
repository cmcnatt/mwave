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
classdef PowerMatrix < IEnergyComp
    % Defines a power matrix for a WEC
    
    properties (Access = private)    
        mat;
        h;
        rho;
        t;
        hs;
        t02;
        specType;
        devCount;
    end
    
    properties (Dependent)
        Matrix;
        Hs;
        T02;
        H;
        T;
        Rho;
        SpecType;
        DeviceCount; 
    end
    
    methods
        
        function [pmat] = PowerMatrix(mat, Hs, T02, T, H, Rho, specType)
            pmat.mat = mat;
            pmat.hs = Hs;
            pmat.t02 = T02;
            pmat.t = T;
            pmat.h = H;
            pmat.rho = Rho;
            pmat.specType = specType;
            pmat.devCount = 1;
        end
        
        function [val] = get.Rho(pmat)
            % fluid density
            val = pmat.rho;
        end
        
        function [val] = get.H(pmat)
            % Get the water depeth
            val = pmat.h;
        end
        
        function [val] = get.T(pmat)
            % wave periods at which the power matrix was computed
            val = pmat.t;
        end
        
        function [val] = get.DeviceCount(pmat)
            % The number of devices evaluated in the EnergyComp
            val = pmat.devCount;
        end
        function [] = set.DeviceCount(pmat, val)
            % The number of devices evaluated in the EnergyComp
            if ~isInt(val)
                error('The DeviceCount must be an integer');
            end
            pmat.devCount = val;
        end
        
        function [val] = get.Matrix(pmat)
            % the array power matrix
            val = pmat.mat;
        end
        
        function [val] = get.Hs(pmat)
            % the signifincant wave heights used
            val = pmat.hs;
        end
        
        function [val] = get.T02(pmat)
            % the T02 used
            val = pmat.t02;
        end
        
        function [val] = get.SpecType(pmat)
            % the type of wave spectrum used to compute the power matrix
            val = pmat.specType;
        end
        
        function [val] = PowerRAO(pmat)
            warning('PowerRAO not defined for PowerMatrix');
            val = [];
        end
        
        function [energy] = AnnualEnergyProd(pmat, waveClim, varargin)
            opts = checkOptions({{'interpPmat'}}, varargin);
            
            % default is to interpolate the wave climate
            interpPmat = opts(1);
            
            if interpPmat
                pmatI = pmat.InterpolateTo(waveClim.Hs('intended'), waveClim.T02('intented'));
                waveClimI = waveClim;
            else
                pmatI = pmat;
                waveClimI = waveClim.InterpolateTo(pmat.hs, pmat.t02);
            end
            
            hrsYr = 24*365;
            freqOccs = waveClimI.FreqOccurance;
            Pow = hrsYr*freqOccs.*pmatI.Matrix; % kWh/yr

            % Total power: MWh/yr
            energy = sum(sum(Pow))./1e3;
        end
        
        function [pow] = AveragePower(pmat, spectrum, varargin)
            % For bretschneider, T02 = 0.71*Tp (T02 =
            % sqrt(m2/m0)*Tp) Ref: Holthuijsen, Waves in
            % Oceanic and Coastal Waters
            hss = spectrum.SigWaveHeight;
            t02s = spectrum.PeakPeriod/0.71;
                
            pmatI = pmat.InterpolateTo(hss, t02s);
            pow = pmatI.Mat;
        end
                
        function [pmatI] = InterpolateTo(pmat, Hs, T)
            [t0M, hsM] = meshgrid(pmat.t02, pmat.hs);
            [TM, HsM] = meshgrid(T, Hs);
            
            mat2 = interp2(t0M, hsM, pmat.mat, TM, HsM);
            
            mat2(isnan(mat2)) = 0;
            
            pmatI = PowerMatrix(mat2, pmat.hs, pmat.t02, ...
                pmat.t, pmat.h, pmat.rho, pmat.specType);
        end
        
        function [] = PlotScatter(pmat, varargin)
            [opts, args] = checkOptions({{'skip', 1}}, varargin);
            
            skip = 2;
            if opts(1)
                skip = args{1};
            end
            
            hs_ = pmat.Hs;
            HsLab = cell(1,length(hs_));
            for n = 1:length(hs_)
                HsLab{n} = num2str(hs_(n), '%4.1f');
            end
            t02_ = pmat.T02;
            T02Lab = cell(1,length(t02_));
            for n = 1:length(t02_)
                T02Lab{n} = num2str(t02_(n), '%4.1f');
            end

            imagesc(pmat.mat);
            indsHs = 1:skip:length(hs_);
            indsT = 1:skip:length(t02_);
            set(gca, 'XAxisLocation','top',...
                'ytick', indsHs, 'yticklabel', HsLab(indsHs),...
                'xtick', indsT, 'xticklabel', T02Lab(indsT));
            fet;
            xlabel('T02 (s)');
            ylabel('Hs (m)');
            cb = colorbar;
            ylabel(cb, 'power (kW)');
        end
    end
    
    methods (Static)
        
        function [pmat] = CreatePowerMatrix(comp, Hs, T, varargin)
            
            if ~isa(comp, 'IEnergyComp')
                error('The comp must be of type IEnergyComp');
            end
            
            [opts, args] = checkOptions({{'waveClim', 1}, {'minPow', 1}, ...
                {'minOcc', 1}, {'specType', 1}, {'makeObj'}}, varargin);
            
            type = 'bretschneider';
            if opts(4)
                type = args{4};
            end
            
            if opts(1)
                waveClim = args{1};
            else
                if ~isempty(comp.H)
                    waveClim = WaveClimate.MakeWaveClimate(type, Hs, T, 1./comp.T, 'H', comp.H, varargin{:});
                else
                    waveClim = WaveClimate.MakeWaveClimate(type, Hs, T, 1./comp.T, varargin{:});
                end
            end
            
            plim = [];
            if opts(2)
                plim = args{2};
            end
            
            occlim = [];
            if opts(3)
                occlim = args{3};
            end
            
            makeObj = opts(5);

            [Mc, Nc] = waveClim.Size;
            freqOccs = waveClim.FreqOccurance;
            Efs = waveClim.EnergyFlux;

            pmat = zeros(Mc, Nc);

            for m = 1:Mc
                for n = 1:Nc
                    if ~isempty(occlim)
                        if (freqOccs(m, n) <= occlim)
                            % ignore sea states that occur less than an hour per year
                            continue;
                        end
                    end

                    if ~isempty(plim)
                        if (sum(Efs{m, n}) <= plim)
                            % ignore sea states below a power threshold
                            continue;
                        end
                    end

                    % power in kW
                    pmat(m, n) = comp.AveragePower(waveClim.WaveSpectra(m, n));
                end
            end
            
            if makeObj
                pmat = PowerMatrix(pmat, waveClim.Hs('intended'), ...
                    waveClim.T02('intended'), ...
                    waveClim.T, waveClim.H, waveClim.Rho, type);
            end
        end
    end
end