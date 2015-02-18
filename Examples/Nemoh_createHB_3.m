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
% In this run, we will create a Flap 'HydroBody' in Nemoh. 
%
% Note - the HydroBody computation using a "CylArray" can only be performed
% with the latest version of Nemoh!
%
% For more info, see 
%   - Nemoh_createHB_1
%   - Wamit_createHB_1
%   - Wamit_createHB_2
%   - hydrBody_1
%   - array_inter_1

%% Set up run

run_name = 'nem_hb_3';         
folder = [mwavePath 'Examples\BemRuns\' run_name];  

rho = 1000;     
h = 50;

len = 5;
beam = 30;
draft = 20;

Nx = 5;            
Ny = 30;                 
Nz = 20;                
wec = FloatingFlap(rho, len, beam, draft, Nx, Ny, Nz);
wec.Modes = ModesOfMotion([0 0 0 0 1 0]);
wec.Handle = 'Flap';

nem_run = NemohRunCondition(folder);    

nem_run.Rho = rho;      
nem_run.H = h;        
nem_run.FloatingBodies = wec;       

Tlims = [8 12];             
nem_run.OmegaLims = 2*pi./fliplr(Tlims);
nem_run.OmegaCount = 5;     


nBetaTot = 40;
nBetaHalf = nBetaTot/2 + 1;
betaLims = [0 pi];
nem_run.BetaLims = betaLims;
nem_run.BetaCount = nBetaHalf;
  
r = 1.1*wec.Rcir;         
nZ = 50;                       
nTheta = 2^6;   

nem_run.CylArray = BemCylArray(r, nTheta, nZ);    
nem_run.ExePath = 'N:\Nemoh\Nemoh_mer\my-nemoh\Release';
nemoh.ComputeHS = false;
nem_run.WriteRun;               % Good to go

%% Run Nemoh
                       
nem_run.Run('Background');           

%% Read results, and save useful objects

nem_result = NemohResult(nem_run);  
nem_result.ReadResult;          

% Even though the field points are set up using a BemCylArray, they are
% still read with WavePoints
waveCir = nem_result.WavePoints;
hydroForces = nem_result.HydroForces;

%% create HydroBody and save it


hydBody = computeHydroBody(waveCir, hydroForces, wec, 'SigFigCutoff', 5,...
    'AccTrim', 'SymX');

save([folder '\nem_hb3_1_hb'], 'hydBody');

%% look at other HydroBody

load 'nem_hb3_1_hb';

