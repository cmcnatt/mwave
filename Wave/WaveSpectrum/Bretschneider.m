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
classdef Bretschneider < WaveSpectrum
    
    methods
        function [ws] = Bretschneider(Hs, Tp, T, varargin)
                % Creates a Bretschneider spectrum. Hs is the significant
                % wave height. Tp is the modal period (s). T is list of periods at which to
                % evaluate the spectrum.  E is the energy density spectrum.
                [opts, args] = checkOptions({{'Cos2s', 3}}, varargin);
                spread = opts(1);

                ws = ws@WaveSpectrum();
                if (nargin > 0)
                    f = 1./T;

                    E = Bretschneider.MakeSpec(Hs, Tp, T, varargin{:});
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
        function [E] = MakeSpec(Hs, Tp, T, varargin)
                % Creates a Bretschneider spectrum. Hs is the significant
                % wave height. Tp is the modal period (s). T is list of periods at which to
                % evaluate the spectrum.  E is the energy density spectrum.
                %
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

                coef = 5/16*fp^4*Hs^2;

                E = coef./(f.^5).*exp(-5/4*(fp./f).^4);
        end
        
        function [To] = ConverterT(Ti, typeIn, typeOut)
            
            % convert to T02
            if strcmpi(typeIn, 'T02')
                T02 = Ti;
            elseif strcmpi(typeIn, 'Te')
                T02 = 1/1.206*Ti;
            elseif strcmpi(typeIn, 'Tp')
                T02 = 0.71*Ti;
            else
                error('typeIn not recognized');
            end
            
            % convert T02 to typeOut.
            if strcmpi(typeOut, 'T02')
                To = T02;
            elseif strcmpi(typeOut, 'Te')
                To = 1.206*T02;
            elseif strcmpi(typeOut, 'Tp')
                To = 1/0.71*T02;
            else
                error('typeOut not recognized');
            end
        end
    end
    
end