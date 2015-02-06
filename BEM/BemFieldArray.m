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
classdef BemFieldArray
    % Defines the field array (gird) of points for a BEM calculation
    
    properties (Access = private)
        start;
        deltas;
        numberPoints;
    end

    properties (Dependent)
        Start;          % The starting location
        Deltas;         % The spacing in each direction
        NumberPoints;   % The number of points in each direction
        Lengths;        % The lengths in each direction
    end
    
    methods (Static, Access = private)
        function [] = checkDeltas(del)
            if (length(del) == 3)    
                for n = 1:3
                    if (del(n) <= 0)
                        error('All deltas must be greater than 0');
                    end
                end
            else
                error('The deltas must be a vector of length 3');
            end
        end
        
        function [np2] = checkNumberPoints(np)
             if (length(np) == 3)   
                np2 = zeros(size(np));
                for n = 1:3
                    if (np(n) <= 0)
                        error('The number of points must be greater than 0');
                    end
                    np2(n) = round(np(n));
                end
            else
                error('The number of points must be a vector of length 3');
            end
        end
    end
    
    methods
        % Constructor 
        function [array] = BemFieldArray(start, deltas, numberPoints)
            % Constructor
            if (nargin == 0)
                array.start = [0 0 0];
                array.deltas = [1 1 1];
                array.numberPoints = [1 1 1];
            else
                if (length(start) == 3)    
                    array.start = start;
                else
                    error('The start location must be a vector of length 3');
                end
                
                array.checkDeltas(deltas);
                array.deltas = deltas;
                
                array.checkNumberPoints(numberPoints);
                array.numberPoints = numberPoints;
            end
        end
        
        % Start
        function [st] = get.Start(array)
            st = array.start;
        end
        function [array] = set.Start(array, st)
            if (length(st) == 3)    
                array.start = st;
            else
                error('The start location must be a vector of length 3');
            end
        end
        
        % Deltas
        function [del] = get.Deltas(array)
            del = array.deltas;
        end
        function [array] = set.Deltas(array, del)
            array.checkDeltas(del);
            array.deltas = del;
        end
        
        % NumberPoints
        function [np] = get.NumberPoints(array)
            np = array.numberPoints;
        end
        function [array] = set.NumberPoints(array, np)
            np2 = array.checkNumberPoints(np);
            array.numberPoints = np2;
        end
        
        function [lens] = get.Lengths(array)
            del = array.deltas;
            np = array.numberPoints;
            
            lens = del.*(np - [1 1 1]);
        end
        
        % Computes meshgrids of the array points
        function[X, Y, Z] = GetArrayPoints(array)
            ends = (array.numberPoints - 1).*array.deltas + array.start;
            x = array.start(1):array.deltas(1):ends(1);
            y = array.start(2):array.deltas(2):ends(2);
            z = array.start(3):array.deltas(3):ends(3);
            
            [X, Y, Z] = meshgrid(x, y, z);
            
            X = squeeze(X);
            Y = squeeze(Y);
            Z = squeeze(Z);
        end
        
        function[points] = GetPointsList(array)
            % returns an Nx3 list of points where N = Nx*Ny*Nz for setting
            % up field arrays with older versions of WAMIT.
            np = array.numberPoints;
            totnp = np(1)*np(2)*np(3);
            points = zeros(totnp, 3);
            
            X0 = array.start;
            del = array.deltas;
            
            pnum = 0;
            
            for l = 1:np(1)
                x = X0(1) + (l-1)*del(1);
                for m = 1:np(2)
                    y = X0(2) + (m-1)*del(2);
                    for n = 1:np(3)
                        z = X0(3) + (n-1)*del(3);
                        
                        pnum = pnum + 1;
                        points(pnum, :) = [x, y, z];
                    end
                end
            end
        end
    end
end