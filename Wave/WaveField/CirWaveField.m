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
classdef CirWaveField < FuncWaveField
    % A wave field provides information (pressure, elevation, velocity,
    % significant wave height, spectrums) at discrete points in space
    % specified on either a grid (array) or at arbitrary (x,y,z) points.  
    %
    % Constructs a wave field of circular (cylindrical) waves radiating 
    % (or scattered) from a source.
    %
    % The inputs are 
    %   - loc: the [x,y] location of the point source
    %   - rho: the water density
    %   - waveComps: the wave components used. i.e. the wave period,
    %   amplitude, phase, water depth, and direction
    %   - coeffs: the coefficients to each wave component.  coeffs is a
    %   cell array where the number of cells is equal to the number of wave
    %   components. Each cell contains a vector or a matrix.  
    %       The coefficients are for the  propagating and evanescent 
    %   circular wave modes.  The matrix should be of size 
    %   (M + 1) x (2N + 1). 
    %       The first row contains the coefficients of the propagating 
    %   modes and each subsequent row contains the coefficients of the 
    %   subsequent evanescent mode. M is the number of evanescent wave 
    %   modes.
    %       The columns are the coefficients for the circular modes.  The 
    %   number of colums is length 2N + 1, where N is the cutoff index in 
    %   the infinite summation  And so the summation goes from -N to N.  
    %   The first  point in the vector is the -N coefficeint. The middle 
    %   point is the  n = 0 coefficient, and the last point is the +N 
    %   coefficient.   
    %   - isarray: indicates whether the wave field is an array or a list
    %   of points
    %   - If isarray, then also must contain a mesh grid of X and Y
    %   - If not is array, then must contain a [Np x 3] vector of {x,y,z}
    %   points
    
    properties (Access = private)
        coeffs;
    end
    
    properties (Dependent)
        Coeffs;         
        KochinFuncs;
    end
    
    methods
        % Constructor
        function [wf] = CirWaveField(rho, waves, isarray, varargin)
            opts = checkOptions({'NoVel', 'Interp'}, varargin);
            noVel = opts(1);
            interp = opts(2);
            
            if (waves.IsPlane)
                error('Input waves must be CirWaves');
            end
            
            if (waves.IsIncident)
                error('Input waves must be outgoing circular waves');
            end
            
            if (interp)
                % for now can only interpolate on z = 0;
                if (isarray)
                    % TODO: ok, because arrays are only X and Y, but if
                    % array extended to include Z, must handle
                else
                    z = varargin{1}(:,3);
                    if (any(z) ~= 0)
                        error('Can only interpolate CirWaveField points for z = 0');
                    end
                end
            end
            
            if (interp)
                if (isarray)
                    args = cell(3,1);
                    args(1:2) = varargin(1:2);
                    args{3} = 'Interp';
                else
                    args = cell(2,1);
                    args(1) = varargin(1);
                    args{2} = 'Interp';
                end
            else
                args = varargin;
            end
            
            wf = wf@FuncWaveField(rho, ~noVel, waves, @CirWaveField.Compute, isarray, args{:});
        end
    end
    
    methods (Static)
        
        function [p, vel] = Compute(rho, compVel, waves, n, x, y, z, varargin)

            loc = waves.Origin;
            if (length(varargin) > 1)
                interp = true;
            else
                interp = false;
            end

            omega = waves.Omega(n);
            h = waves.H;
            k = waves.K(n);
            g = waves.G;

            x = x - loc(1);
            y = y - loc(2);
            
            A = waves.Coefs{n};
            
            [L, M, A0, Am, Anm] = CirWaves.SortCoefs(A);
            
            if (L > 0)
                kl = IWaves.SolveForK(omega, h, g, 'Evanescent', L);
                kl = kl(2:end);
            end
            
            expimtheta = cell(1,M);
            expNimtheta = cell(1,M);
            
            if (interp)
                mx = max(max(abs(x)));
                my = max(max(abs(y)));
                R0 = sqrt(mx^2 + my^2);
                dx = min(min(x(2:end, 2:end) - x(1:end-1, 1:end-1)));
                dy = min(min(y(2:end, 2:end) - y(1:end-1, 1:end-1)));

                dr = min([dx dy]);
                r = (0:dr:R0)';

                dtheta = 2*asin(dy/R0);
                theta = -pi-dtheta:dtheta:pi+dtheta;

                Nr = length(r);
                Ntheta = length(theta);
                
                z = zeros(Nr, Ntheta);
            else
                r = sqrt(x.^2 + y.^2);
                theta = atan2(y,x);   
            end
            
            kr = k*r;                
            
            % Propagating modes
            H0 = besselh(0,2,kr);
            H1 = besselh(1,2,kr);
                
            if (interp)
                p = A0(1)*H0*ones(1, Ntheta);
                if (compVel)
                    ur = A0(1)*H1*ones(1, Ntheta);
                    utheta = zeros(Nr, Ntheta);
                end
            else
                p = A0(1).*H0;
                if (compVel)
                    ur = A0(1).*H1;
                    utheta = 0;
                end
            end

            % higher components
            Hmm1 = H0;
            Hm = H1;

            for m = 1:M
                Hmp1 = besselh(m+1,2,kr);

                cmtheta = cos(m*theta);
                smtheta = sin(m*theta);
                expimtheta{m} = cmtheta + 1i*smtheta;
                expNimtheta{m} = cmtheta - 1i*smtheta;

                Tm = Am(1,m)*expimtheta{m};
                Tnm = Anm(1,m)*expNimtheta{m};

                if (interp)
                    p = p + Hm*(Tm + (-1)^m*Tnm);
                
                    if (compVel)
                        ur = ur + 0.5*(Hmm1 - Hmp1)*(Tm + (-1)^m*Tnm);
                        utheta = utheta + m*Hm*(Tm - (-1)^m*Tnm);

                        Hmm1 = Hm;
                    end
                else
                    p = p + Hm.*(Tm + (-1)^m*Tnm);

                    if (compVel)
                        ur = ur + 0.5*(Hmm1 - Hmp1).*(Tm + (-1)^m*Tnm);
                        utheta = utheta + m*Hm.*(Tm - (-1)^m*Tnm);

                        Hmm1 = Hm;
                    end
                end

                Hm = Hmp1;
            end
            
            if (interp)
                % for interp - z is assumed zero
                p = rho*g*p;
                
                if (compVel)
                    w = 1i*g*k/omega*tanh(k*h).*p;
                    ur = -1i*g*k/omega.*ur;
                    utheta = -g/omega.*utheta;
                end
            else
                if (isinf(h))
                    Kp = exp(k*z);

                    Kpw = Kp;
                else
                    Kp = cosh(k*(h + z))./cosh(k*h);

                    if (compVel)
                        Kpw = sinh(k*(h + z))./cosh(k*h);
                    end
                end

                p = rho*g*Kp.*p;

                if (compVel)
                    w = 1i*g*k/omega*Kpw.*p;
                    ur = -1i*g*k/omega*Kp.*ur;
                    utheta = -g/omega*Kp.*utheta;
                end
            end
            
            % evanescent modes 
            % TODO: still need add velocity computation for evanescent modes
            for l = 1:L
                K0 = besselk(0,kl(l)*r);
                
                coskl = cos(kl(l)*(z + h));
                
                if (interp)
                    p = p + rho*g*coskl.*A0(l+1).*K0*ones(1,Ntheta);
                else
                    p = p + rho*g*coskl.*A0(l+1).*K0;
                end
                
                for m = 1:M
                    Km = besselk(m,kl(l)*r);
                    
                    if (interp)
                        p = p + rho*g*coskl.*Km*(Am(l+1,m)*expimtheta{m} + Anm(l+1,m)*expNimtheta{m});
                    else
                        p = p + rho*g*coskl.*Km.*(Am(l+1,m)*expimtheta{m} + Anm(l+1,m)*expNimtheta{m});
                    end
                end
            end
            
            if (interp)
                [Theta, R] = meshgrid(theta, r);
                ri = sqrt(x.^2 + y.^2);
                thetai = atan2(y,x);  

                p = interp2(R.', Theta.', p.', ri, thetai);
                if (compVel)
                    ur = interp2(R.', Theta.', ur.', ri, thetai);
                    utheta = interp2(R.', Theta.', utheta.', ri, thetai);
                    w = interp2(R.', Theta.', w.', ri, thetai);
                    
                    theta = thetai;
                end
            end
            
            if (compVel)
                cth = cos(theta);
                sth = sin(theta);
                u = ur.*cth - utheta.*sth;
                v = ur.*sth + utheta.*cth;

                vel = zeros([3, size(x)]);

                switch (ndims(vel))
                    case 2
                        vel(1,:) = u;
                        vel(2,:) = v;
                        vel(3,:) = w;
                    case 3
                        vel(1,:,:) = u;
                        vel(2,:,:) = v;
                        vel(3,:,:) = w;
                    case 4
                        vel(1,:,:,:) = u;
                        vel(2,:,:,:) = v;
                        vel(3,:,:,:) = w;
                    otherwise
                        error('velociy in plane wave has a unknown dimensions');
                end
            else
                vel = [];
            end           
        end
    end
end