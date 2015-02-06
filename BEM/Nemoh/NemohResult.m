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
classdef NemohResult < IBemResult
    % Reads Wamit output files
    %
    % It can be set up by using the NemohRunConditions used to create the 
    % run as the constructor argument.  
    
    % Otherwise, it requires the path to the folder in which the output 
    % files are located and it 
    
    methods

        function [result] = NemohResult(runCondition)
            % Constructor
            if (nargin == 0)
                result.folder = ' ';
                result.runName = [];
            else 
                if (isa(runCondition, 'NemohRunCondition'))
                    result.rho = runCondition.Rho;
                    result.folder = runCondition.Folder;
                    result.floatingbodies = runCondition.FloatingBodies;
                    result.dof = HydroBodyComp.GetDoF(result.floatingbodies);
                    result.t = runCondition.T;
                    result.nT = length(result.t);
                    result.beta = runCondition.Beta;
                    result.nB = length(result.beta);
                    result.h = runCondition.H;
                    %result.fieldPoints = runCondition.FieldPoints;
                    result.fieldArray = runCondition.FieldArray;
                    result.fieldPoints = [];
                    result.cylArray = runCondition.cylArray;
                    if (~isempty(result.fieldPoints) || ~isempty(result.fieldArray))
                        result.solveField = true;
                    end
                else
                    error('Constructor input must be of type NemohRunCondition')
                end
            end
            result.hasBeenRead = 0;
        end

        function [] = ReadResult(result, varargin)
            % Reads the Wamit results

            useSing = false;
            removeBodies = false;
            for l = 1:length(varargin)
                if (strcmp(varargin{l}, 'UseSingle'))
                    useSing = true;
                end
                if (strcmp(varargin{l}, 'RemoveBodies'))
                    removeBodies = true;
                end
            end
            
            % Count
            nbod = length(result.floatingbodies);
            nT = result.nT;
                                    
            % Gravity, depth and density        
            result.g = 9.80660; % TODO: read this value
                                     
            % Excitation forces
            % Need to fix for multiple directions
            nB = result.nB;
            fid = fopen([result.folder,'\Results\ExcitationForce.tec'],'r');
            fgetl(fid);             % read header lines
            dofs = zeros(1, nbod);
            for m = 1:nbod
                dofs(m) = result.floatingbodies(m).Modes.DoF;
                for l = 1:dofs(m)
                    fgetl(fid);     % read header lines
                end
            end

            doft = sum(dofs);
            Famp = zeros(nT, nB, doft);
            Fphi = zeros(nT, nB, doft);

            for l = 1:nB
                fgetl(fid);             % read header lines
                for m = 1:nT
                    ligne = fscanf(fid,'%f', 1 + 2*doft);
                    % omega(n) = ligne(1);
                    for n = 1:doft
                        Famp(m, l, n) = ligne(2*n);
                        Fphi(m, l, n) = ligne(2*n+1);
                    end;
                end;
                fgetl(fid);         % for some reason need to read a blank line here
            end
             
            % Note that the complex conjuage is taken here because Nemoh
            % computes for a time depedence of exp(-i*omega*t), while
            % mwave assumes a time depedent of exp(i*omega*t);
            Fe = Famp.*exp(-1i*Fphi);
            fclose(fid);
             
            % Radiation forces
            fid = fopen([result.folder,'\Results\RadiationCoefficients.tec'],'r');
            fgetl(fid);             % read header lines
            for m = 1:nbod
                for l = 1:dofs(m)
                    fgetl(fid);     % read header lines
                end
            end

            A = zeros(nT, doft, doft);
            B = zeros(nT, doft, doft);
            
            for l = 1:doft
                fgetl(fid);         % read header lines
                for m = 1:nT
                    ligne = fscanf(fid,'%f',1 + 2*doft);
                    for n = 1:doft
                        A(m, l, n) = ligne(2*n);
                        B(m, l, n) = ligne(2*n + 1);
                    end;
                    fgetl(fid);
                end;
            end;
            fclose(fid);
            
            % Hydrostatic matrix
            C = zeros(doft, doft);
            
            fid = fopen([result.folder,'\Mesh\KH.dat'],'r');
            
            if (fid > 0)
                for n = 1:doft   
                    ligne = fscanf(fid,'%g %g', doft);
                    C(n,:)=ligne;
                end;            
                fclose(fid);
            end
            
            result.hydroForces = HydroForces(result.t, result.beta, A, B, C, Fe, result.h, result.rho);
                        
            % Wave array
            if (~isempty(result.fieldArray))
                nProb = nT*(doft + nB);
                fsfiles = dir([result.folder,'\Results\freesurface.*.dat']);
                fsfiles = {fsfiles.name};
                
                if (length(fsfiles) ~= nProb)
                    error('The number of free surface files does not match the expected number of problems');
                end
                
                probn = 0;
                
                rho = result.rho;
                g = result.g;
                h = result.h;
                t = result.t;
                beta = result.beta;
                
                EtaS = cell(nT, nB);
                EtaR = cell(nT, doft);
                
                
                % not sure about the body order
                for l = 1:nbod
                    for m = 1:nT
                        % Diffraction
                        for n = 1:nB
                            probn = probn + 1;
                            [eta, xy] = nemoh_readFieldPoints([result.folder,'\Results\' fsfiles{probn}]);
                            %eta = conj(eta);
                            [Eta, X, Y] = result.reshapeGrid(eta, xy);
                            
                            EtaS{m, n} = Eta;
                        end
                        
                        % Radiation
                        for n = 1:dofs(l)
                            probn = probn + 1;
                            [eta, xy] = nemoh_readFieldPoints([result.folder,'\Results\' fsfiles{probn}]);
                            %eta = conj(eta);
                            [Eta, X, Y] = result.reshapeGrid(eta, xy);
                            
                            EtaR{m, n} = Eta;
                        end
                    end
                end
                
                rwfs(result.dof, 1) = WaveField;
                V = [];
                
                for m = 1:doft
                    P = zeros([nT, size(X)]);
                    for n = 1:nT
                        P(n,:,:) = rho*g*EtaR{n, m};
                    end
                    
                    rwfs(m) = WaveField(rho, g, h, t, P, V, 1, X, Y);
                end
                rwavefield = WaveFieldCollection(rwfs, 'MotionIndex', (1:result.dof));
                
                swfs(nB, 1) = WaveField;
                iwfs(nB, 1) = WaveField;
                
                for m = 1:nB
                    P = zeros([nT, size(X)]);
                    for n = 1:nT
                        P(n,:,:) = rho*g*EtaS{n, m};
                    end
                    
                    swf = WaveField(rho, g, h, t, P, V, 1, X, Y);
                    wcs = PlaneWaves(ones(size(t)), t, beta(m)*ones(size(t)), h);
                    iwf = PlaneWaveField(rho, wcs, 1, X, Y);
                    iwfs(m) = iwf;
                    swfs(m) = swf;
                end
                
                iwavefield = WaveFieldCollection(iwfs, 'Direction', beta);
                swavefield = WaveFieldCollection(swfs, 'Direction', beta);

                result.waveArray = FBWaveField(iwavefield, swavefield, rwavefield);
                    
            end
            % Wave points
            
            result.hasBeenRead = 1;
        end
    end
    
    methods (Access = private)
        
        function [Eta, X, Y] = reshapeGrid(result, eta, xy)
            nx = length(unique(xy(:,1)));
            ny = length(unique(xy(:,2)));
            X = reshape(xy(:,1), nx, ny);
            Y = reshape(xy(:,2), nx, ny);
            Eta = reshape(eta, nx, ny);
        end
        
    end

end