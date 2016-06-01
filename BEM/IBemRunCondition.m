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
classdef IBemRunCondition < handle
    % BemRunCondition interface and abstract class.  
    % Creates input files to a BEM (WAMIT or Nemoh) run
    
    properties (SetAccess = protected, GetAccess = protected)        
        rho;
        h;
        exePath;
        folder;
        floatBods;
        geoFiles;
        fieldArray;
        cylArray;
        incZeroFreq;
        incInfFreq;
    end
    
    properties (Dependent)
        Rho;                % Fluid density 
        H;                  % Water depth   
        ExePath;            % Path location of BEM executables
        Folder;             % Folder location of run (input/output files) 
        FieldArray;         % Array of field points (Input is a BemFieldArray)
        CylArray;           % Cylindrical array of field poinst (Input is a BemCylArray)
        IncZeroFreq;        % Include the zero frequency added mass solution
        IncInfFreq;         % Include the infinite freqeuncy added mass solution
    end
    
    properties (Abstract, Dependent)
        T;                  % The wave periods (get only)
        Beta;               % The wave directions (get only)
        FloatingBodies;     % FloatingBodies (plural because of multiple bodies) (Nemoh?)
    end
    
    methods (Abstract)
        WriteRun(run, varargin) % Writes the input files for a BEM run
        Run(run, varargin)      % Runs the executable
    end
    
    methods (Abstract, Access = protected)   
        makeFolderAndSub(run, fold);
    end
    
    methods 
        
        function [rh] = get.Rho(run)
            % Get the water density
            rh = run.rho;
        end
        function [run] = set.Rho(run, rh)
            % Set the water density
            if (rh > 0)
                run.rho = rh;
            else
                error('Rho must be a positive number');
            end
        end
        
        function [h_] = get.H(run)
            % Get the water depeth
            h_ = run.h;
        end
        function [run] = set.H(run, h_)
            % Set the water depth
            if (h_ > 0)
                run.h = h_;
            else
                error('The depth must be a positive number');
            end
        end
        
        function [ep] = get.ExePath(run)
            % Path location of BEM executables
            ep = run.exePath;
        end
        function [run] = set.ExePath(run, ep)
            % Path location of BEM executables
            if (ischar(ep))
                run.exePath = ep;
            else
                error('ExePath must be a string');
            end
        end
        
        function [fol] = get.Folder(run)
            % Get the folder location where the input and output files will go
            fol = run.folder;
        end
        function [run] = set.Folder(run, fol)
            % Set the folder location where the input and output files will go
            if (ischar(fol))
                run.folder = fol;
                run.makeFolderAndSub(fold);
            else
                error('Folder must be a string');
            end
        end
        
        function [fa] = get.FieldArray(run)
            % Get the field point array to be evaluated
            fa = run.fieldArray;
        end
        function [run] = set.FieldArray(run, fa)
            % Set the field point array to be evaluated
            if(~isa(fa, 'BemFieldArray'))
                error('Field array must be of type BemFieldArray');
            end
            run.fieldArray = fa;
        end
        
        function [ca] = get.CylArray(run)
            % Get the cylinder point array to be evaluated
            ca = run.cylArray;
        end
        function [run] = set.CylArray(run, ca)
            % Set the cylinder point array to be evaluated
            if(~isa(ca, 'BemCylArray'))
                error('CylArray array must be of type BemCylArray');
            end
            run.cylArray = ca;
        end
        
        function [val] = get.IncZeroFreq(run)
            % Include the zero frequency added mass solution
            val = run.incZeroFreq;
        end
        function [] = set.IncZeroFreq(run, val)
            % Include the zero frequency added mass solution
            if (~isBool(val))                
                error('IncZeroFreq must be a boolean');
            end
            run.incZeroFreq = val;
        end
        
        function [val] = get.IncInfFreq(run)
            % Include the infinite frequency added mass solution
            val = run.incInfFreq;
        end
        function [] = set.IncInfFreq(run, val)
            % Include the infinite frequency added mass solution
            if (~isBool(val))                
                error('IncInfFreq must be a boolean');
            end
            run.incInfFreq = val;
        end
    end
end