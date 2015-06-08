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
classdef CirWaveComp 
    
    methods (Static)
        function [am] = PlaneWaveCoef(M, beta)
            am = exp(-1i*(-M:M)*(beta+pi/2));
        end
        function [am, A, theta] = CoefsFromSpread(M, s, varargin)
            % Produces the curved (or short-crested) wave coefficient for curved
            % incident waves from a cos squared spreading defined by s: 
            %
            %   A^2 = cos(1/2*theta).^(2*s)
            %   
            %   where A is the amplitude density spectrum.. 
            % 
            % Inputs:   M - coefficient limit (number of coefs is 2*M+1)
            %           s - spreading factor 
            %           k - the wave number
            %           'RandPhase' - optional input, to create a random phase
            %
            % Ouputs:   am - circular wave coefficients
            %           A - the directional spread
            %           theta - theta points for B

            opts = checkOptions({'Rand'}, varargin);
            rand = opts(1);

            % directional coefficients
            Nthet = 4096;
            dthet = 2*pi/Nthet;
            %theta = 0:dthet:(2*pi-dthet);
            theta = -pi:dthet:(pi-dthet);

            A2 = cos(1/2*theta).^(2*s);
            Amp = trapz(theta,A2);
            A2 = A2./Amp;
            
            A = sqrt(A2);
            Amp = trapz(theta,A2);
            A = A./Amp;

            if (rand)
                [C, L] = CirWaveComp.CovMatrixFromFunc(A2, 2*M);
                x = randn(2*M+1,1);
                am = L*x;
                am = am.';
            else
                am = CirWaveComp.CoefsFromFunc(A, M);
            end
        end
        
        function [C, L] = CovMatrixFromFunc(A, M)
%             s = 40;
%             M = 50;
%             
%             dbeta = pi/2000;
%             beta = -pi:dbeta:(pi - dbeta);
%             Nbeta = length(beta);
%             
%             S = cos(1/2*beta).^(2*s);
%             intS = trapz(beta, S);
%             
%             S = S./intS;
%             fori = pi/Nbeta*fft(S);
%             
%             Cm = ones(2*M+1,1);
%             
%             Cm(M+1) = fori(1);
%             Cm(M+2:2*M+1) = fori(2:M+1);
%             Cm(M:-1:1) = fori(Nbeta:-1:Nbeta-M+1);
%             
%             rndval = 10^-10;
%             
%             Cm = rndval*round(real(Cm)./rndval) + 1i*rndval*round(imag(Cm)./rndval);
%             
%             mm = (-M:M).';
%             
%             Cm = (-1i).^mm.*Cm;

            Cm = CirWaveComp.CoefsFromFunc(A, M);
            
            N = floor(M/2);
            C = zeros(2*N+1);
                    
            for m = -N:N
                for n = -N:N
                    ind = m-n;
                    C(N+m+1, N+n+1) = Cm(M+ind+1);
                end
            end
            
            [L, p] = chol(C, 'lower');
        end
        
        function [am] = CoefsFromFunc(A, M)
            N = length(A);
            
            fori = 2*pi/N*fft(A);
            
            am = ones(2*M+1,1);
            am(M+1) = fori(1);
            am(M+2:2*M+1) = fori(2:M+1);
            am(M:-1:1) = fori(N:-1:N-M+1);

            m = (-M:M).';
            % here +1i is used instead of -1i because the input is given
            % from -pi to pi, rather than 0 to 2*pi, which is what the fft
            % expects
            am = exp(1i*m*pi/2).*am;
            am = am.';
        end
        
        function [MatOut] = Resize2M(MatIn, M)
            [nrow, ncol] = size(MatIn);
            
            Min = ncol;
            if (mod(Min,2) ~= 1)
                error('The number of cols of the input matrix must be odd');
            end
            Min = (Min - 1)/2;
            
            if (nrow == 1)
                if (M > Min)
                    MatOut = zeros(1, 2*M+1);
                    MatOut((M + 1 - Min):(M + 1 + Min)) = MatIn;
                elseif (M < Min)
                    MatOut = MatIn((Min + 1 - M):(Min + 1 + M));
                else
                    MatOut = MatIn;
                end
            elseif (nrow == ncol)
                if (M > Min)
                    MatOut = zeros(2*M+1,2*M+1);
                    MatOut((M + 1 - Min):(M + 1 + Min),(M + 1 - Min):(M + 1 + Min)) = MatIn;
                elseif (M < Min)
                    MatOut = MatIn((Min + 1 - M):(Min + 1 + M), (Min + 1 - M):(Min + 1 + M));
                else
                    MatOut = MatIn;
                end
            else
                error('Input matrix not of recognizable size');
            end
        end
    end
    
end