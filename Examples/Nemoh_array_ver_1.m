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

%% Load the Nemoh HydroBodies

run_name = 'nem_hb_1';     
load([mwavePath 'Examples\BemRuns\' run_name '\nem_hb1_1_hb']);

cylHydBod = hydBody;

run_name = 'nem_hb_2';    
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

tic
Aa = arrayComp.A;        
acomptime = toc;
fprintf('Array computation time: %5.3f s\n', acomptime);
Ba = arrayComp.B;        
Ca = arrayComp.C;        
Fexa = arrayComp.Fex;    

x = -100:2:100;
y = x;
[X, Y] = meshgrid(x, y);

wfa = arrayComp.WaveField(true, X, Y, 'NoVel');  

%% Load Wamit results for same case

load wam_array_ver_res

%% Or rerun Wamit

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

wam_run.Run;           

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

save([folder '\wam_array_ver_res'], 'wfw', 'wcomp', 'Aw', 'Bw', 'Cw', 'Fexw');

%% Compare results

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
title('Array comp (Nemoh)');

subplot(2,1,2);
pcolor(X,Y,abs(etaw{iT}));
fet;
set(gca, 'clim', clims);
title('Wamit');

%%
figure;
i1 = 1;
i2 = 7;
subplot(5,2,1);
plot(T, [squeeze(Aa(:,i1,i2)) squeeze(Aw(:,i1,i2))]);
title('Surge c1-surge c2');
subplot(5,2,2);
plot(T, [squeeze(Ba(:,i1,i2)) squeeze(Bw(:,i1,i2))]);

i1 = 2;
i2 = 8;
subplot(5,2,3);
plot(T, [squeeze(Aa(:,i1,i2)) squeeze(Aw(:,i1,i2))]);
title('Sway c1-sway c2');
subplot(5,2,4);
plot(T, [squeeze(Ba(:,i1,i2)) squeeze(Bw(:,i1,i2))]);

i1 = 3;
i2 = 9;
subplot(5,2,5);
plot(T, [squeeze(Aa(:,i1,i2)) squeeze(Aw(:,i1,i2))]);
title('Heave c1-heave c2');
subplot(5,2,6);
plot(T, [squeeze(Ba(:,i1,i2)) squeeze(Bw(:,i1,i2))]);

i1 = 2;
i2 = 13;
subplot(5,2,7);
plot(T, [squeeze(Aa(:,i1,i2)) squeeze(Aw(:,i1,i2))]);
title('Surge c1-Flap');
subplot(5,2,8);
plot(T, [squeeze(Ba(:,i1,i2)) squeeze(Bw(:,i1,i2))]);

i1 = 7;
i2 = 13;
subplot(5,2,9);
plot(T, [squeeze(Aa(:,i1,i2)) squeeze(Aw(:,i1,i2))]);
title('Surge c2-Flap');
subplot(5,2,10);
plot(T, [squeeze(Ba(:,i1,i2)) squeeze(Bw(:,i1,i2))]);

%%
figure;
i1 = 1;
subplot(5,1,1);
plot(T, abs([squeeze(Fexa(:,1,i1)) squeeze(Fexw(:,1,i1))]));
title('Surge c1');

i1 = 7;
subplot(5,1,2);
plot(T, abs([squeeze(Fexa(:,1,i1)) squeeze(Fexw(:,1,i1))]));
title('Surge c2');

i1 = 3;
subplot(5,1,3);
plot(T, abs([squeeze(Fexa(:,1,i1)) squeeze(Fexw(:,1,i1))]));
title('Heave c2');

i1 = 11;
subplot(5,1,4);
plot(T, abs([squeeze(Fexa(:,1,i1)) squeeze(Fexw(:,1,i1))]));
title('Pitch c2');

i1 = 13;
subplot(5,1,5);
plot(T, abs([squeeze(Fexa(:,1,i1)) squeeze(Fexw(:,1,i1))]));
title('Flap');







