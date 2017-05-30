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
% In this run, we will create a cylinder 'HydroBody' in Nemoh. 
%
% Note - the HydroBody computation using a "CylArray" can only be performed
% with the latest version of Nemoh!
%
% For more info on HydroBodies, see 
%   - Wamit_createHB_1
%   - Wamit_createHB_2
%   - hydrBody_1
%   - array_inter_1

%% Set up run

run_name = 'nem_hb_1';         
folder = [mwavePath 'Examples\BemRuns\' run_name];  

rho = 1000;     
h = 50;

% create a cylinder - In mwave, NemohRunConditions takes the same 
% FloatingBody inputs WamitRunCondition does. However, Nemoh is not 
% currently set up to do Generalized modes.  
rad = 10;
hei = 15;
draft = 10;

Ntheta = 48;            
Nr = 8;                 
Nz = 24;                
wec = FloatingCylinder(rho, rad, hei, draft, Ntheta, Nr, Nz);


nem_run = NemohRunCondition(folder);    % The NemohRunCondition just takes 
                                        % a folder that is the exact
                                        % location where the run goes

nem_run.Rho = rho;      
nem_run.H = h;        
nem_run.FloatingBodies = wec;       

% The NemohRunCondition creates a Nemoh.cal file, which is used by the
% preProcessor. The Nemoh.cal file takes radial frequency limits and count
Tlims = [8 12];             % 8 to 12 seconds
nem_run.OmegaLims = 2*pi./fliplr(Tlims);
nem_run.OmegaCount = 5;     % 5 frequencies - evenly spaced in omega

% The wave directions are the same as the frequency
nBeta = 40;
betaLims = [0 2*pi*(1-1./nBeta)];
nem_run.BetaLims = betaLims;
nem_run.BetaCount = nBeta;
  
                                        
% In Wamit_createHB_1, the cylindrical array of points was created
% manually. However, NemohRunCondition (and Nemoh) cannot take a list of 
% points. The cylindrical array of points is specified in Nemoh.cal as one 
% line: 
%
%   r   nTheta  nZ
%
% Note - this is only available in the lastest version of Nemoh
%
% This is setup in NemohRunCondition using a BemCylArray

r = 1.1*wec.Rcir;                                                                            
nZ = 50;                       
nTheta = 2^6;   

nem_run.CylArray = BemCylArray(r, nTheta, nZ);    
nem_run.ExePath = 'N:\Nemoh\Nemoh_mer\my-nemoh\Release';
nem_run.WriteRun;               % Good to go

%% Run Nemoh

nem_run.Run('Background');           

%% Read results, and save useful objects

nem_result = NemohResult(nem_run);  
nem_result.ReadResult;          

% Even though the field points are set up using a BemCylArray, they are
% still read with WavePoints
waveCir = nem_result.WavePoints;
hydroForces = nem_result.FreqDomForces;

%% create HydroBody and save it


hydBody = computeHydroBody(waveCir, hydroForces, wec, 'SigFigCutoff', 5,...
    'AccTrim');

save([mwavePath 'Examples\HydroBodies\nem_hb1_1_hb'], 'hydBody');

