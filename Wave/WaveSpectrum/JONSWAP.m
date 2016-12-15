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
classdef JONSWAP < WaveSpectrum
    
    methods
        function [ws] = JONSWAP(Hs, Tp, T, gamma, varargin)
                % Creates a JONSWAP spectrum. Hs is the significant
                % wave height. Tp is the modal period (s). T is list of 
                % periods at which to evaluate the spectrum.  gamma is the 
                % spectral spreading parameter.
                [opts, args] = checkOptions({{'Cos2s', 3}}, varargin);
                spread = opts(1);

                ws = ws@WaveSpectrum();
                if (nargin > 0)
                    f = 1./T;

                    E = JONSWAP.MakeSpec(Hs, Tp, T, gamma, varargin{:});
                    if (spread)
                        s = args{1}{1};
                        dirc = args{1}{2};
                        dir = args{1}{3};

                        G = cosSpectSpread(s, dirc, dir);
                        E = E.'*G;
                        ws.setSpectrum(E, f, dir);
                    else
                        ws.setSpectrum(E, f);
                    end
                end
        end
    end
    
    methods (Static)
        function [E] = MakeSpec(Hs, Tp, T, gamma, varargin)
            % Creates a JONSWAP spectrum. Hs is the significant
            % wave height. Tp is the modal period (s). T is list of 
            % periods at which to evaluate the spectrum.  gamma is the 
            % spectral spreading parameter. 
            % E is the variance density spectrum (m^2/Hz).

            [opts] = checkOptions({'T02'}, varargin);
            isT02 = opts(1);

            if (isT02)
                % For bretschneider, T02 = 0.71*Tp (T02 =
                % sqrt(m2/m0)*Tp) Ref: Holthuijsen, Waves in
                % Oceanic and Coastal Waters
                fp = 1./(Tp/0.71);
            else
                fp = 1./Tp;
            end
            f = 1./T;

            fcomp = 0.001:0.001:2;
            HsParamLims = [0 Hs];
            E = JONSWAP.makeSpec(Hs, fp, fcomp, gamma);
            HsOut = 4*sqrt(trapz(fcomp, E));
            func = Hs - HsOut;
            funcLims = [Hs func];

            N = 100;
            err = 0.01;

            for n = 1:N
                HsParam = mean(HsParamLims);
                E = JONSWAP.makeSpec(HsParam, fp, fcomp, gamma);
                HsOut = 4*sqrt(trapz(fcomp, E));
                func = Hs - HsOut;
                if abs(func) < err
                    break;
                end

                if func < funcLims(1) && func > 0
                    HsParamLims(1) = HsParam;
                    funcLims(1) = func;
                elseif func > funcLims(2)
                    HsParamLims(2) = HsParam;
                    funcLims(2) = func;
                else
                    error('Not converging');
                end
            end
            
            E = JONSWAP.makeSpec(HsParam, fp, f, gamma);
        end
    end
    
    methods (Access = private, Static)
        function [E] = makeSpec(Hs, fp, f, gamma)
            % http://www.wikiwaves.org/Ocean-Wave_Spectra
            coef = 5/16*fp^4*Hs^2;
            
            sigma = [0.07 0.09];
            
            inds{1} = f <= fp;
            inds{2} = f > fp;

            E = zeros(size(f));
            
            for n = 1:2
                E(inds{n}) = coef./(f(inds{n}).^5).*exp(-5/4*(fp./f(inds{n})).^4)...
                    .*gamma.^exp(-(f(inds{n})-fp).^2./(2*sigma(n)*fp).^2);
            end
        end
    end
    
end