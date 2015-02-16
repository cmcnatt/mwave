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
% In this example, we'll use two HydroBodies created in Nemoh to perform
% and array computation, and then compare it to the array computed directly
% with Nemoh.
%
% The HydroBodies were computed with
%   - Nemoh_createHB_1
%   - Nemoh_createHB_2
%
% For more info see:
%   - array_inter_1
%   - array_inter_2

%% Load the HydroBodies

run_name = 'nem_hb_1';  % Or we can be more explicit if necessary         
load([mwavePath 'Examples\BemRuns\' run_name '\nem_hb1_1_hb']);

cylHydBod = hydBody;

run_name = 'nem_hb_2';  % Or we can be more explicit if necessary         
load([mwavePath 'Examples\BemRuns\' run_name '\nem_hb2_1_hb']);

flapHydBod = hydBody;

%% Compute the array interaction with mwave

% Cylinders (C) and flap (F) at the following positions
%
%   C(-30,40)
%           
%               F(30, 0)
%
%   C(-30,-40)



% Create an array of HydroBodies
hbs(1) =  HydroBody(cylHydBod);    
hbs(2) =  HydroBody(cylHydBod); 
hbs(3) = HydroBody(flapHydBod);

% The positions
hbs(1).XYpos = [-30, 40];
hbs(2).XYpos = [-30, -40];
hbs(3).XYpos = [30, 0];

% Initialize the array computation with unit amplitude plane waves
T = hydBody.T;         
h = hydBody.H;      
a = ones(size(T)); 
beta = 0;

iwaves = PlaneWaves(a, T, beta, h);

% create the initial array computation object
arrayComp = HydroArrayComp(hbs, iwaves);

Aa = arrayComp.A;        
Ba = arrayComp.B;        
Ca = arrayComp.C;        
Fexa = arrayComp.Fex;    

x = -200:2:200;
y = x;
[X, Y] = meshgrid(x, y);

wfa = arrayComp.WaveField(true, X, Y, 'NoVel');  

%% Compute the same array with Nemoh

run_name = 'nem_array_1';         
folder = [mwavePath 'Examples\BemRuns\' run_name];  

nem_run = NemohRunCondition(folder);    

nem_run.Rho = hbs(1).Rho;      
nem_run.H = hbs(1).H;        
nem_run.FloatingBodies = hbs;       

T = hbs.T;        
nem_run.OmegaLims = 2*pi./[T(1) T(end)];
nem_run.OmegaCount = length(T);
nem_run.BetaLims = [beta beta];
nem_run.BetaCount = 1;

nps = [length(x) length(y) 1];
lens = [(x(end) - x(1)) (y(end) - y(1)) 1];
dels = lens./(nps - [1 1 0]);
starts = -lens./2;
starts(3) = 0;

nem_run.FieldArray = BemFieldArray(starts, dels, nps);
nem_run.ExePath = 'N:\Nemoh\Nemoh_mer\my-nemoh\Release';
nem_run.WriteRun;               
                        
nem_run.Run('Background');           

%% Read Nemoh results

nem_result = NemohResult(nem_run);  
nem_result.ReadResult;          

% Even though the field points are set up using a BemCylArray, they are
% still read with WavePoints
wfn = nem_result.WaveArray;
hydroForces = nem_result.HydroForces;

An = hydroForces.A;        
Bn = hydroForces.B;        
Cn = hydroForces.C;        
Fexn = hydroForces.Fex;   

%% And to double check, compute with Wamit

run_name = 'wam_array_1';         
folder = [mwavePath 'Examples\BemRuns\' run_name];  

wam_run = WamitRunCondition(folder, run_name);    

wam_run.Rho = hbs(1).Rho;      
wam_run.H = hbs(1).H;        
wam_run.FloatingBodies = hbs;       

wam_run.T = hbs(1).T;
wam_run.Beta = 0;

nps = [length(x) length(y) 1];
lens = [(x(end) - x(1)) (y(end) - y(1)) 1];
dels = lens./(nps - [1 1 0]);
starts = -lens./2;
starts(3) = 0;

wam_run.FieldArray = BemFieldArray(starts, dels, nps);
wam_run.WriteRun;               
                        
wam_run.Run('Background');           

%% Read Wamit results

wam_result = WamitResult(wam_run);  
wam_result.ReadResult;          

% Even though the field points are set up using a BemCylArray, they are
% still read with WavePoints
wfw = wam_result.WaveArray;
hydroForces = wam_result.HydroForces;
wcomp = HydroBodyComp(hydroForces, wam_run.FloatingBodies);

Aw = hydroForces.A;        
Bw = hydroForces.B;        
Cw = hydroForces.C;        
Fexw = hydroForces.Fex;   

save([folder '\wam_array_ver_res'], 'wfw', 'wcomp');

%%

type = 'Diffracted';
iT = 4;
clims = [0.8, 1.2];

etaa = wfa.Elevation(type);
etaw = wfw.Elevation(type);

figure;
subplot(2,1,1);
pcolor(X,Y,abs(etaa{iT}));
fet;
set(gca, 'clim', clims);

subplot(2,1,2);
pcolor(X,Y,abs(etaw{iT}));
fet;
set(gca, 'clim', clims);






