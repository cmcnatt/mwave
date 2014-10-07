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
classdef CirWaves < IWaves
    
    properties (Access = private)
        origin;
        mlim;
        el;
        isInc;
        kochinFs;
    end
    
    properties (Dependent)
        IsIncident;      % Indicates whether it is an incident wave.
        IsPlane;         % Indicates whether it is a long-crested plane wave  
        Mlim;            % Truncation values of the circular coefficient for the circular modes. 
        L;               % Truncation values of the circular coefficients for the evanescent modes. 
        KochinFunc;      % The Kochin function is the far-field amplitude as a function of direction.
        Coefs;          % Circular wave field coefficients
        Origin;         % Location of the origin of the circular waves
    end
    
    methods
        
        function [wav] = CirWaves(inOrOut, org, coefs, t, h)
            % Constructor
            wav.nT = -1;
            if (nargin ~= 0)
                a = ones(size(t));
                beta = zeros(size(t));
                wav.init(a, t, beta, h);
            end
            
            if (strcmp(inOrOut, 'In'))
                wav.isInc = true;
            elseif (strcmp(inOrOut, 'Out'))
                wav.isInc = false;
            else
                error('First input to contstructor must be either ''In'' (for incident waves) or ''Out'' (for radiated and scattered waves)');
            end
            
            if all(size(org) ~= [1 2])
                error('Origin must be a 1x2 array');
            end
            
            if (~iscell(coefs))
                error('The coefficients must be supplied as a cell array');
            end
            
            if any(size(coefs) ~= [wav.nT 1])
                error('The coefficients must be of size NT x 1');
            end
 
            mli = zeros(wav.nT, 1);
            l = zeros(wav.nT, 1);
            for m = 1:wav.nT;
                [L, M, ~, ~, ~] = CirWaves.SortCoefs(coefs{m});
                mli(m) = M;
                l(m) = L;
            end
            
            if (wav.isInc && any(l ~= 0))
                error('An incident circular wave cannot have evanescent components');
            end
            
            wav.mlim = mli;
            wav.el = l;
            wav.coefs = coefs;
            wav.origin = org;
        end
        
        function [ii] = get.IsIncident(wav)
            % Indicates whether it is an incident wave.
            ii = wav.isInc;
        end
        
        function [ip] =  get.IsPlane(wav)
            % Indicates whether it is a long-crested plane wave  
            ip = false;
        end
        
        function [mli] = get.Mlim(wav)
            % Truncation values of the circular coefficient for the circular modes. 
            mli = wav.mlim;
        end
                  
        function [l] = get.L(wav)
            % Truncation values of the circular coefficients for the evanescent modes. 
            l = wav.el;
        end
        
        function [kfs] = get.KochinFunc(wav)
            % The Kochin function is the far-field amplitude as a function of direction.
            % Gets the Kochin Functions associated with this wave field.
            if (isempty(wav.kochinFs))
                kfs(wav.nT, 1) = KochinFunc;
                for n = 1:wav.nT
                    kfs(n) = KochinFunc('Coefs', wav.coefs{n});
                end
                wav.kochinFs = kfs;
            else
                kfs = wav.kochinFs;
            end
                
        end
        
        function [A] = get.Coefs(wav)
            % Circular wave field coefficients
            A = wav.coefs;
        end
        
        function [or] = get.Origin(wav)
            % Location of the origin of the circular waves
            or = wav.origin;
        end
        
        function [as] = IncAmps(wav, M, varargin)
            % Circular coefficients describing the wave.
            
            if (wav.isInc)
                func = @(x,y) (besselj(x, y));
            else
                func = @(x,y) (besselj(x, y) - 1i*bessely(x, y));
            end
            
            isorg = true;
            iwav = 0;
            if (~isempty(varargin))
                pos = varargin{1};
                if (any(abs(wav.origin - pos) >= 1e-9*[1, 1]))
                    isorg = false;
                end
                
                if (length(varargin) > 1)
                    iwav = varargin{2};
                end
            else
                pos = wav.origin;
            end
            
            if (iwav == 0)
                startm = 1;
                stopm = wav.nT;
                as = cell(wav.nT, 1);
            else
                startm = iwav;
                stopm = startm;
            end

            for m = startm:stopm
                a = zeros(1, 2*M+1);

                Mf = wav.mlim(m);

                thisA = wav.coefs{m};

                if (Mf <= M)
                    a(M-Mf+1:M+Mf+1) = thisA(1,:);
                else
                    a = thisA(1,Mf-M+1:Mf+M+1);
                end

                if (~isorg)
                    Tij = CirWaves.basisT_(func, M, pos, wav.origin, wav.k(m));
                    a = Tij.'*(a.');
                    a = a.';
                end
                
                a = wav.a(m)*exp(1i*wav.epsilon(m))*a;

                if (iwav == 0)
                    as{m} = a;
                else
                    as = a;
                end
            end
        end
    end
    
    methods (Static)
        
        function [L, M, A0, Am, Anm] = SortCoefs(A)
            [nrow, ncol] = size(A);

            if (rem(ncol,2) ~= 1)
                error('There must be an odd number of circular coefficients (rows)');
            end
            
            M = (ncol - 1)/2;
            L = nrow - 1;

            % negative indexed coefficients come first but are ordered from
            % -infty to -1, and so need to be flipped.
            Anm = fliplr(A(:,1:M));
            
            A0 = A(:,M+1);
            Am = A(:,M+2:(2*M+1));
        end
        
        function [Tij] = BasisTH(M, Ci, Cj, k0)
            func = @(x,y) (besselj(x, y) - 1i*bessely(x, y));
            Tij = CirWaves.basisT_(func, M, Ci, Cj, k0);
        end
         
        function [Tij] = BasisTJ(M, Ci, Cj, k0)
            func = @(x,y) (besselj(x, y));
            Tij = CirWaves.basisT_(func, M, Ci, Cj, k0);
         end
    end
           
    methods (Static, Access = private)
        
        function [Tij] = basisT_(func, M, Ci, Cj, k0)
             
             L = sqrt((Ci(1) - Cj(1))^2 + (Ci(2) - Cj(2))^2);
             k0L = k0*L;
             
             alpha = atan2(Ci(2) - Cj(2), Ci(1) - Cj(1));
             
             Tij = zeros(2*M+1, 2*M+1);
                          
             for m = -M:M
                 for n = -M:M
                     %H = besselh(m-n, 2, k0L);
                     F = func(m-n, k0L); 
                     
                     Tij(m+M+1, n+M+1) = F*exp(1i*(m-n)*alpha);
%                      if (abs(H - H1)/abs(H1) > 1e-7)
%                          error('hankel function not the same');
%                      end
                 end
             end
         end
    end
end