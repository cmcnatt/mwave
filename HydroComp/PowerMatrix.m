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
        f;
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
        Se;
        T02;
        H;
        T;
        Rho;
        SpecType;
        DeviceCount;
    end
    
    methods
        
        function [pmat] = PowerMatrix(mat, Hs, T02, T, H, Rho, specType)
            % Note: it is not straightforward to add an 'Se' option because
            % of the use of nargin here, so we cannot use varargin with
            % checkOpts as we might in other codes. 
            %   If wish to make
            % a power matrix based on an Se vector, the user should compute
            % the Hs matrix that corresponds to the Se vector, and input
            % that to this constructor. (i.e. the reverse of the code on
            % get.Se)

            pmat.f =[];
            pmat.mat = mat;
            pmat.hs = Hs;
            pmat.t02 = T02;
            if nargin > 3
                pmat.t = T;
                pmat.h = H;
                pmat.rho = Rho;
                pmat.specType = specType;
            else
                pmat.t = [];
                pmat.h = [];
                pmat.rho = [];
                pmat.specType = [];
            end
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
            % the significant wave heights used
            val = pmat.hs;
        end
        
        function [val] = get.Se(pmat)
            % the wave steepnesses used
            Hs = pmat.hs;
            H = Hs*2*sqrt(2)/4; % Find wave height of that regular wave that has equivalent energy to the irregular wave of Hs.
            
            % Find Te from T02 values stored
            TtypeIn = 'T02'; TtypeOut = 'Te';
            Te = Bretschneider.ConverterT(pmat.t02, TtypeIn, TtypeOut);
            
            if size(Hs,1) > 1 && size(Hs,2) > 1
                L = repmat(1.5608*Te.*Te,size(Hs,1),1); % Find wavelength of regular wave with period Te
                SeAll = H./L; % Find wave steepnesses
                Se = SeAll(:,1); % Pick out just one column - all should be identical if Hs matrix was generated using a single Se vector.
            else
                Harray = repmat(H',1,length(Te));
                Larray = repmat(1.5608*Te.*Te,length(Hs),1);
                Se = Harray./Larray;
            end
            
            val = Se';
        end
        
        function [val] = get.T02(pmat)
            % the T02 used
            val = pmat.t02;
        end
        
        function [val] = get.SpecType(pmat)
            % the type of wave spectrum used to compute the power matrix
            val = pmat.specType;
        end
        
        function [p] = PowerAt(pmat, Hs, T, varargin)
            
            [opts, args] = checkOptions({{'T02'}, {'Tp'}, {'Te'}}, varargin);
            isT02 = opts(1);
            isTp = opts(2);
            isTe = opts(3);
            
            if isempty(pmat.f)
                pmat.f =  griddedInterpolant({pmat.hs, pmat.t02}, pmat.Matrix, 'linear', 'none');
            end
            
            Ti = T;
            if isT02
                Ti = T;
            elseif isTp
                Ti = Bretschneider.ConvertT(T, 'Tp', 'T02');
            elseif isTe
                Ti = Bretschneider.ConvertT(T, 'Te', 'T02');
            end
            
            p = pmat.f(Hs, Ti);
            p(isnan(p)) = 0;
        end
        
        function [val] = PowerRAO(pmat)
            warning('PowerRAO not defined for PowerMatrix');
            val = [];
        end
        
        function [energy] = AnnualEnergyProd(pmat, waveClim, varargin)
            opts = checkOptions({{'interpPmat'}}, varargin);
            
            % default is to interpolate the wave climate
            interpPmat = opts(1);
            
            % Check power matrix and wave climate objects are of the same type.
            if ( sum(size(pmat.Hs)>1) == 2 && sum(size(waveClim.Hs('intended'))>1) ~= 2 ) ...
                    || ( sum(size(pmat.Hs)>1) ~= 2 && sum(size(waveClim.Hs('intended'))>1) == 2 )
                error('In order to compute AEP, power matrix object and wave climate object must be of same type - either Hs-T or Se-T.')
            end
            
            typeSe = 0;
            if sum(size(pmat.Hs)>1) == 2 % i.e. if using steepness instead of Hs to parameterise the power matrix
                typeSe = 1;
            end
            
            if interpPmat
                if typeSe
                    pmatI = pmat.InterpolateTo(waveClim.Hs('intended'), waveClim.T02('intended'), 'typeSe', waveClim.Se('Intended'));
                else
                    pmatI = pmat.InterpolateTo(waveClim.Hs('intended'), waveClim.T02('intended'));
                end
                waveClimI = waveClim;
            else
                pmatI = pmat;
                if typeSe
                    waveClimI = waveClim.InterpolateTo(pmat.hs, pmat.t02, 'typeSe', pmat.Se);
                else
                    waveClimI = waveClim.InterpolateTo(pmat.hs, pmat.t02);
                end
            end
            
            hrsYr = 24*365;
            freqOccs = waveClimI.FreqOccurance;
            Pow = hrsYr*freqOccs.*pmatI.Matrix; % kWh/yr
            
            % Total power: MWh/yr
            energy = sum(sum(Pow))./1e3;
        end
        
        function [pow] = AveragePower(pmat, varargin)
            % For bretschneider, T02 = 0.71*Tp (T02 =
            % sqrt(m2/m0)*Tp) Ref: Holthuijsen, Waves in
            % Oceanic and Coastal Waters
            
            if nargin == 2
                spectrum = varargin{1};
                hss = spectrum.SigWaveHeight;
                ti = spectrum.PeakPeriod;
                ttype = 'Tp';
            elseif nargin == 3
                hss = varargin{1};
                ti = varargin{2};
                ttype = 'T02';
            elseif nargin == 4
                hss = varargin{1};
                ti = varargin{2};
                ttype = varargin{3};
            end
            
            pow = pmat.PowerAt(hss, ti, ttype);
            
            %             pmatI = pmat.InterpolateTo(hss, t02s);
            %             pow = pmatI.Mat;
        end
        
        function [pmatI] = InterpolateTo(pmat, Hs, T, varargin)
            [opts, args] = checkOptions({{'typeSe',1}}, varargin);
            
            if opts(1)
                [t0M, hsM] = meshgrid(pmat.t02, pmat.Se);
                HsMatrix = Hs;
                Se = args{1};
                [TM, HsM] = meshgrid(T, Se);
            else
                [t0M, hsM] = meshgrid(pmat.t02, pmat.hs);
                [TM, HsM] = meshgrid(T, Hs);
            end
            
            mat2 = interp2(t0M, hsM, pmat.mat, TM, HsM);
            
            mat2(isnan(mat2)) = 0;
            
            pmatI = PowerMatrix(mat2, pmat.hs, pmat.t02, ...
                pmat.t, pmat.h, pmat.rho, pmat.specType);
        end
        
        function [] = PlotScatter(pmat, varargin)
            [opts, args] = checkOptions({{'skip', 1},{'compareWith',1}}, varargin);
            
            skip = 2;
            if opts(1)
                skip = args{1};
            end
            if opts(2) % i.e. wish to make comparison plot
                pmat2 = args{2};
            end
            
            indsHs = 1:skip:length(pmat.Hs);
            indsT = 1:skip:length(pmat.T02);
            
            if opts(2)
                plotScatter(pmat.T02, pmat.Hs, 100*pmat.mat./pmat2.mat, 'xinds', indsT, 'yinds', indsHs);
            else
                plotScatter(pmat.T02, pmat.Hs, pmat.mat, 'xinds', indsT, 'yinds', indsHs);
            end
            
            xlabel('T02 (s)');
            ylabel('Hs (m)');
            cb = colorbar;
            if opts(2)
                ylabel(cb, 'power (% difference)');
            else
                ylabel(cb, 'power (kW)');
            end
        end
    end
    
    methods (Static)
        
        function [pmat, idptos, errs, Dptos, Dpars, tdas, pmatMech, pmatAC] = CreatePowerMatrix(comp, Hs, T, varargin)
            
            if ~isa(comp, 'IEnergyComp')
                error('The comp must be of type IEnergyComp');
            end
            
            isSpec = isa(comp, 'SpecDomComp');
            isTime = isa(comp, 'TimeDomComp');
            
            [opts, args] = checkOptions({{'waveClim', 1}, {'minPow', 1}, ...
                {'minOcc', 1}, {'specType', 1}, {'makeObj'}, ...
                {'ratedPow', 1}, {'dptos', 1}, {'HsLim', 1}, {'seed',1}, {'ParallelOn'}, ...
                {'fullMatrix'}, {'powMech'}, {'powAC'}, {'beta', 1}}, varargin);
            
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
            
            ratedPow = [];
            if opts(6)
                ratedPow = args{6};
            end
            
            dptos = [];
            if opts(7)
                dptos = args{7};
            end
            hslim = [];
            if opts(8)
                hslim = args{8};
            end
            
            if opts(9)
                seed = args{9};
            else
                seed = []; % Will just use random number generator to choose a seed if one isn't selected by the user
            end
            
            if isTime && opts(10)
                tdParOn = 1; % Choose best parallel programming option with time domain computation - i.e. use parsim instead of parfor loop.
            else
                tdParOn = 0;
            end
            
            fullMatrix = 0;
            if opts(11)
                fullMatrix = 1;
            end
            
            mechPowMat = 0;
            if opts(12)
                mechPowMat = 1;
                if ~isTime
                    error('Mechanical power matrix can only be generated using a time-domain model.')
                end
            end
            acPowMat = 0;
            if opts(13)
                acPowMat = 1;
                if ~isTime
                    error('AC power matrix can only be generated using a time-domain model.')
                end
            end

            if opts(14)
                beta0 = args{14}; % Allow input of a non-zero angle between WEC and waves
            else
                beta0 = 0; % otherwise, assume WEC is aligned perfectly to waves
            end
            
            [Mc, Nc] = waveClim.Size;
            
            Te = Bretschneider.ConverterT(waveClim.T02('intended'), 't02', 'te');
            Hs = waveClim.Hs('intended');
            if size(Hs,1) > 1 && size(Hs,2) > 1
                typeSe = 1; % when wave steepness is used, there is a unique Hs for each Se-Te pairing
            else
                typeSe = 0;
            end
            
            if ~isempty(hslim)
                if ~typeSe
                    Mc = indexOf(Hs, hslim); % If using Se, Hs does not simply trim some rows off the matrix, so keep Mc x Nc in size.
                end
            end
            
            freqOccs = waveClim.FreqOccurance;
            Efs = waveClim.EnergyFlux;
            
            pmat = zeros(Mc, Nc);
            pmatMech = zeros(Mc, Nc);
            pmatAC = zeros(Mc, Nc);
            idptos = ones(Mc, Nc);
            errs = zeros(Mc, Nc);
            Dptos = cell(Mc, Nc);
            Dpars = cell(Mc, Nc);
            tdas = cell(Mc, Nc);
            
            if isSpec
                Dpto0 = comp.FreqDomComp.Dpto;
            else
                Dpto0 = comp.Dpto;
            end
            
            if tdParOn % if (TD with parsim)
                
                if fullMatrix
                    
                    for n = 1:Nc
                        % power in kW
                        if ~isempty(dptos)
                            error('When ParallelOn option is used, code is not yet setup to use multiple dptos.');
                        else
                            comp.SetDpto(Dpto0);
                            
                            tic
                            [pmat(:, n), tdasTemp, pmatMech(:,n), pmatAC(:,n)] = comp.AveragePower(waveClim.WaveSpectra(:, n),'seed',seed,'ParallelOn','beta',beta0);
                            for m = 1:length(tdasTemp)
                                tdas{m,n} = tdasTemp{m,1};
                            end
                            
                            fprintf('\nn = %i/%i, Te = %4.1f, run time = %4.1f s\n', ...
                                n, Nc, Te(n), toc);
                        end
                        
                        if ~isempty(dptos)
                            Dptos{m, n} = dptos{ind};
                        end
                        
                    end
                    
                else
                    
                    if ~isempty(dptos)
                        error('When ParallelOn option is used, code is not yet setup to use multiple dptos.');
                    else
                        comp.SetDpto(Dpto0);
                        
                        % To ensure CPUs are maximally utilised, convert
                        % matrices to vectors. Also enforce Hs limit and
                        % occurrence threshold.
                        if isempty(hslim)
                            hslim = Inf;
                        end
                        if ~typeSe
                            HsMatrix = repmat(Hs',1,Nc);
                        else
                            HsMatrix = Hs;
                        end
                        runIndsMatrix = (waveClim.freqOcc.*(HsMatrix <= hslim) > occlim);
                        runInds = (waveClim.freqOcc(HsMatrix <= hslim) > occlim);
                        seaStateSpectra = waveClim.WaveSpectra(runInds);
                        
                        tic
                        [pmatVector, tdasTemp, pmatMechVector, pmatACVector] = comp.AveragePower(seaStateSpectra,'seed',seed,'ParallelOn','beta',beta0);
                        
                        % Insert vector of tda objects back into matrix
                        % form.
                        tdaCounter = 1;
                        for j = 1:size(runIndsMatrix,2) % Must loop down the columns first as if reading a book top to bottom along the lines.
                            for i = 1:size(runIndsMatrix,1)
                                if runIndsMatrix(i,j) == 1
                                    tdas{i,j} = tdasTemp{tdaCounter,1};
                                    tdaCounter = tdaCounter + 1;
                                end
                            end
                        end
                        
                        % Insert vector of non-zero power matrix entries back into
                        % the power matrix itself.
                        pmat(runInds) = pmatVector;
                        if mechPowMat
                            pmatMech(runInds) = pmatMechVector;
                        end
                        if acPowMat
                            pmatAC(runInds) = pmatACVector;
                        end
                        
                        fprintf('\nTotal run time = %4.1f s\n', toc);
                    end
                    
                end
                
            else % if (FD) or (SD) or (TD with parfor)
                
                parfor m = 1:Mc % this line can be parfor/for
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
                        
                        if typeSe
                            if (Hs(m, n) >= hslim)
                                % ignore sea states for which Hs is greater
                                % than the HsLim
                                continue;
                            end
                        end
                        
                        % power in kW
                        if ~isempty(dptos)
                            powmn = zeros(length(dptos),1);
                            powACmn = zeros(length(dptos),1);
                            powMechmn = zeros(length(dptos),1);
                            
                            errmn = zeros(length(dptos),1);
                            tdamn = cell(length(dptos),1);
                            for o = 1:length(dptos)
                                comp.SetDpto(dptos{o});
                                
                                tic
                                if isSpec
                                    [powmn(o), errmn(o)] = comp.AveragePower(waveClim.WaveSpectra(m, n),'seed',seed);
                                elseif isTime
                                    error('Simulink model may not initialise correctly using parfor - instead, it should be parallelised using parsim.');
                                    [powmn(o), tdamn{o}, powACmn(o), powMechmn(o)] = comp.AveragePower(waveClim.WaveSpectra(m, n),'seed',seed,'beta',beta0);
                                else
                                    [powmn(o)] = comp.AveragePower(waveClim.WaveSpectra(m, n),'seed',seed);
                                end
                                
                                if typeSe
                                    fprintf('\nm = %i/%i, n = %i/%i, o = %i/%i, Hs = %4.1f, Te = %4.1f, run time = %4.1f s\n', ...
                                        m, Mc, n, Nc, o, length(dptos), Hs(m,n), Te(n), toc);
                                else
                                    fprintf('\nm = %i/%i, n = %i/%i, o = %i/%i, Hs = %4.1f, Te = %4.1f, run time = %4.1f s\n', ...
                                        m, Mc, n, Nc, o, length(dptos), Hs(m), Te(n), toc);
                                end
                            end
                            [pmat(m, n), ind] = max(powmn);
                            if acPowMat
                                [pmatAC(m, n), ind] = max(powACmn);
                            end
                            if mechPowMat
                                [pmatMech(m, n), ind] = max(powMechmn);
                            end
                            idptos(m, n) = ind;
                            errs(m, n) = errmn(ind);
                            tdas(m, n) = tdamn(ind);
                            
                            fprintf('\nm = %i/%i, n = %i/%i, Best Dpto Ind = %i\n', ...
                                m, Mc, n, Nc, ind);
                        else
                            comp.SetDpto(Dpto0);
                            tic
                            if isSpec
                                [pmat(m, n), errs(m, n)] = comp.AveragePower(waveClim.WaveSpectra(m, n),'seed',seed);
                            elseif isTime
                                [pmat(m, n), tdas{m, n}, pmatAC(m, n), pmatMech(m, n)] = comp.AveragePower(waveClim.WaveSpectra(m, n),'seed',seed,'beta',beta0);
                            else
                                [pmat(m, n)] = comp.AveragePower(waveClim.WaveSpectra(m, n),'seed',seed);
                            end
                            
                            if typeSe
                                fprintf('\nm = %i/%i, n = %i/%i, Hs = %4.1f, Te = %4.1f, run time = %4.1f s\n', ...
                                    m, Mc, n, Nc, Hs(m,n), Te(n), toc);
                            else
                                fprintf('\nm = %i/%i, n = %i/%i, Hs = %4.1f, Te = %4.1f, run time = %4.1f s\n', ...
                                    m, Mc, n, Nc, Hs(m), Te(n), toc);
                            end
                        end
                        
                        if isSpec
                            Dptos{m, n} = comp.FreqDomComp.Dpto;
                            Dpars{m, n} = comp.FreqDomComp.Dpar;
                        elseif isTime
                            if ~isempty(dptos)
                                Dptos{m, n} = dptos{ind};
                            end
                        else
                            if ~isempty(dptos)
                                Dptos{m, n} = dptos{ind};
                                Dpars{m, n} = comp.Dpar;
                            end
                        end
                        
                        if ~isempty(ratedPow)
                            pmat(m, n) = min([pmat(m, n) ratedPow]);
                            if acPowMat
                                pmatAC(m, n) = min([pmatAC(m, n) ratedPow]);
                            end
                            if mechPowMat
                                pmatMech(m, n) = min([pmatMech(m, n) ratedPow]);
                            end
                        end
                    end
                end
                
            end
            
            % for a spectral domain computation
            % set power for sea-states that failed to converge to 0
            if isSpec
                pmat(errs > 10^-4) = 0;
            end
            
            if makeObj
                pmat = PowerMatrix(pmat, waveClim.Hs('intended'), ...
                    waveClim.T02('intended'), ...
                    waveClim.T, waveClim.H, waveClim.Rho, type);
                if mechPowMat
                    pmatMech = PowerMatrix(pmatMech, waveClim.Hs('intended'), ...
                    waveClim.T02('intended'), ...
                    waveClim.T, waveClim.H, waveClim.Rho, type);
                else
                    pmatMech = [];
                end
                if acPowMat
                    pmatAC = PowerMatrix(pmatAC, waveClim.Hs('intended'), ...
                    waveClim.T02('intended'), ...
                    waveClim.T, waveClim.H, waveClim.Rho, type);
                else
                    pmatAC = [];
                end
            end
        end
    end
end