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
classdef KochinFunc
    
    properties (Access = private)
        isCoefs;
        coefs;
        thet;
        f;
        df;
    end
    
    methods
        function [kf] = KochinFunc(varargin)
            
            if (isempty(varargin))
                kf.isCoefs = [];
                kf.coefs = [];
                kf.thet = [];
                kf.f = [];
                kf.df = [];
            else
                [opts, args] = checkOptions({{'Elevation', 4}, {'Coefs', 1}}, varargin);

                if ((opts(1) && opts(2)) || (~opts(1) && ~opts(2)))
                    error('Must be either an Elevation or Coefs type Kochin function');
                end

                if (opts(1))
                    kf.isCoefs = false;
                    kf.coefs = [];

                    k = args{1}{1};
                    r = args{1}{2};
                    theta = args{1}{3};
                    eta = args{1}{4};

                    kr = k*r;
                    F = sqrt(kr)*exp(1i*kr)*eta;

                    % check to see if theta is equally spaced
                    dtheta = diff(theta);
                    eqlSp = true;
                    for n = 2:length(dtheta)
                        if (abs(dtheta(n) - dtheta(1)) > 1e-12)
                            eqlSp = false;
                            break;
                        end
                    end

                    N = length(F);
                    f2 = zeros(1, N+2);
                    f2(1) = F(N);
                    f2(2:N+1) = F;
                    f2(N+2) = F(1);

                    theta2 = zeros(1, N+2);
                    theta2(1) = theta(N) - 2*pi;
                    theta2(2:N+1) = theta;
                    theta2(N+2) = theta(1) + 2*pi;

                    if (eqlSp)
                        dtheta = dtheta(1);

                        dF = (f2(3:N+2) - f2(1:N))./(2*dtheta);
                    else
                        dF = zeros(1, N);

                        for n = 2:N+1
                            x = theta2(n);
                            x0 = theta2(n-1);
                            x1 = theta2(n);
                            x2 = theta2(n+1);

                            dF(n-1) = f2(n-1)*(2*x - x1 - x2)/((x0 - x1)*(x0 - x2)) + f2(n)*(2*x - x0 - x2)/((x1 - x0)*(x1 - x2)) + f2(n+1)*(2*x - x0 - x1)/((x2 - x0)*(x2 - x1));
                        end
                    end

                    df2 = zeros(1, N+2);
                    df2(1) = dF(N);
                    df2(2:N+1) = dF;
                    df2(N+1) = dF(1);


                    kf.thet = theta2;
                    kf.f = f2;
                    kf.df = df2;
                else
                    kf.isCoefs = true;
                    kf.coefs = args{2};
                    kf.thet = [];
                    kf.f = [];
                    kf.df = [];
                end
            end
        end
        
        function [f_, df_] = Evaluate(kf, theta)
            if (kf.isCoefs)
                A = kf.coefs;
                [~, M, A0, Am, Anm] = CirWaves.SortCoefs(A);

                % just get the propogating coefs
                A0 = A0(1);
                Am = Am(1,:);
                Anm = Anm(1,:);

                f_ = A0*ones(size(theta));
                df_ = zeros(size(theta));

                for m = 1:M
                    Amexp = Am(m)*exp(1i*m*(theta+pi/2));         
                    Anmexp = Anm(m)*exp(-1i*m*(theta+pi/2));
%                      Amexp = Am(m)*exp(1i*m*(theta-pi/2));         % This is the version that corresponds to the wave spectrum
%                      Anmexp = Anm(m)*exp(-1i*m*(theta-pi/2));

                    f_ = f_ + (Amexp + Anmexp);
                    df_ = df_ + m*(1i^(m + 1)*Amexp - 1i^(-m + 1)*Anmexp);
                end

                f_ = 1/pi*f_;                                        % This is the normalization used in the thesis
                df_ = 1/pi*df_;
                
%                 f_ = sqrt(2/pi)*exp(1i*pi/4)*f_;
%                 df_ = sqrt(2/pi)*exp(1i*pi/4)*df_;
            else
                f_ = interp1(kf.thet, kf.f, theta);
                df_ = interp1(kf.thet, kf.df, theta);
            end
        end
    end
end
