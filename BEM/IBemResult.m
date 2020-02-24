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
classdef IBemResult < handle
    % Hold BEM output files
    %
    % It can be set up by using an IBemRunConditions used to create the 
    % run as the constructor argument.  
    
    % Otherwise, it requires the path to the folder in which the output 
    % files are located and it 
    
    properties (Access = protected)
        rho;
        g;
        h;
        folder;
        floatingbodies;
        t;
        nT;
        beta;
        nB;
        dof;
        fieldPoints;
        fieldArray;
        cylArray;
        solveField;
        hydroForces;
        wavePoints;
        waveArray;
        hasBeenRead;
        errLog;
        driftForces; % This drift force structure could contain more than just the mean drift forces that will be computed initially
    end

    properties (Dependent)
        Rho;                % Fluid density
        G;                  % Gravitational constant
        H;                  % Depth
        Folder;             % Folder location of run (string)
        FloatingBodies;     % FloatingBodies (plural because of multiple bodies)
        T;                  % Periods (s)
        Beta;               % Directions (deg)
        FreqDomForces;        % Hydrodynamic forces computed by WAMIT
        WavePoints;         % Wave field at field points
        WaveArray;          % Wave field at array locations
        ErrorLog;
        DriftForces; 
    end
    
    methods (Abstract)
        ReadResult(result, varargin)    % Reads the output files for a BEM run
    end

    methods
        
        function [fol] = get.Folder(result)
            % Get the folder location where the output files are found
            fol = result.folder;
        end
        function [] = set.Folder(result, fol)
            % Set the folder location where the output files are found
            if (ischar(fol))
                result.folder = fol;
            else
                error('Folder must be a string');
            end
            result.hasBeenRead = 0;
        end
        
        function [fbs] = get.FloatingBodies(result)
            % Get the array of floating bodies
            fbs = result.floatingbodies;
        end
        function [] = set.FloatingBodies(result, fbs)
            % Set the array of floating bodies
            for n = 1:length(fbs)
                if (~isa(fbs(n), 'FloatingBody'))
                    error('All Floating Bodies must be of type FloatingBody');
                end
            end
            
            result.floatingbodies = fbs;
            result.hasBeenRead = 0;
        end
        
        function [rh] = get.Rho(result)
            % The water density
            if (result.hasBeenRead)
                rh = result.rho;
            else
                rh = [];
            end
        end
        
        function [g_] = get.G(result)
            % The gravitational constant
            if (result.hasBeenRead)
                g_ = result.g;
            else
                g_ = [];
            end
        end
        
        function [t_] = get.T(result)
            % The wave periods
            if (result.hasBeenRead)
                t_ = result.t;
            else
                t_ = [];
            end
        end

        function [bet] = get.Beta(result)
            % The wave headings
            if (result.hasBeenRead)
                bet = result.beta;
            else
                bet = [];
            end
        end
               
        function [h_] = get.H(result)
            % The water depth
            if (result.hasBeenRead)
                h_ = result.h;
            else
                h_ = [];
            end
        end
        
        function [hf] = get.FreqDomForces(result)
            % The FreqDomForces object
            if (result.hasBeenRead)
                hf = result.hydroForces;
            else
                hf = [];
            end
        end
        
        function [df] = get.DriftForces(result)
            % The DriftForces structure
            if (result.hasBeenRead)
                df = result.driftForces;
            else
                df = [];
            end
        end
        
        function [wp] = get.WavePoints(result)
            % Wave field points, not in an array
            if (result.hasBeenRead)
                wp = result.wavePoints;
            else
                wp = [];
            end
        end
        
        function [wa] = get.WaveArray(result)
            % Wave field points in an array
            if (result.hasBeenRead)
                wa = result.waveArray;
            else
                wa = [];
            end
        end
        
        function [er] = get.ErrorLog(result)
            if (result.hasBeenRead)
                er = result.errLog;
            else
                er = [];
            end
        end
        
    end
    
end