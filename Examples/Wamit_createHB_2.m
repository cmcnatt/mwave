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
% In this run, we will create a 'HydroBody' in Wamit. The HydroBody will be
% a cylinder in heave. For more details on creating a HydroBody, see 
% example Wamit_createHB_1.
%
% The HydroBody created here is used in example hydroBody_2

%% Set up run

run_name = 'wam_hb_2';         
folder = [mwavePath '\Examples\BemRuns\' run_name];  

rho = 1025;     % water density
h = 8;          % water depth
T = 1:0.5:8;    % wave periods (8 s wave is 64 m long at 8 m depth, if 
                % this were scaled to full scale by mulitplying by 5 (40 m
                % water depth), the wave would be 320 m long, and have a
                % period of 17 s)

% set up cylinder dimension
rad = 1;
draft = 1;
hei = 2;

Ntheta = 24;   
Nr = 10;
Nz = 10;

wec = FloatingCylinder(rho, rad, hei, draft, Ntheta, Nr, Nz);
wec.Modes = ModesOfMotion([0 0 1 0 0 0]);   % compute for heave only

N = 40;                         % 40 directions
Beta = 0:2*pi/N:2*pi*(1-1/N);   % from 0 to 2pi

% In Wamit_createHB_1, the cylindrical array of points was created
% manually. However, in WamitRunCondition, one may also use BemCylArray, 
% which automatically creates the points with a cos spacing in z. For
% NemohRunCondition, only BemCylArray is valid (i.e. it doesn't take a list
% of points) - see Nemoh_createHB_1

r = wec.Rcir;                   
dr = 0.1;                       
r = r + dr;                     

nZ = 200;                       
nTheta = 2^8;  

cylArray = BemCylArray(r, nTheta, nZ);                 

wam_run = WamitRunCondition(folder, run_name);  

wam_run.Rho = rho;      
wam_run.T = T;                
wam_run.Beta = Beta;              
wam_run.H = h;        
wam_run.FloatingBodies = wec;       

wam_run.CylArray = cylArray;   
wam_run.WriteRun;               

%% Run Wamit

% wam_run.Run;                           
wam_run.Run('Background');           

%% Read results, and save useful objects

wam_result = WamitResult(wam_run);  
wam_result.ReadResult;          

waveCir = wam_result.WavePoints;
hydroForces = wam_result.HydroForces;


%% create HydroBody and save it

hydBody = computeHydroBody(waveCir, hydroForces, wec, 'SigFigCutoff', 5,...
    'AccTrim');

save([mwavePath 'Examples\HydroBodies\wam_hb2_1_hb'], 'hydBody');
% this HydroBody is used in example hydroBody_2

