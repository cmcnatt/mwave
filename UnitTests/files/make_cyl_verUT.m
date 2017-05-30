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
%% 1 - Compute wave field for hydrobody

% Write WAMIT run
a = 1;
lam = 10;

N = 40;

Beta = 0:2*pi/N:2*pi*(1-1/N);

dr = 0.1;
nZ = 200;
nTheta = 2^7;

mainPath = fileparts(which('mwave.home'));
mainPath = [mainPath '\UnitTests\files\'];

% Circ wave fields - write runs

% Cylinder
name = 'hb';
path = [mainPath name];

h = 10;     
rho = 1000;

rad = 1;
hei = 1.5;
draft = 1;

Ntheta = 64;
Nr = 8;
Nz = 12;
    
wec1 = FloatingCylinder(rho, rad, hei, draft, Ntheta, Nr, Nz);

wec1.Handle = 'cyl';     
wec1.Modes = ModesOfMotion([1 1 1 1 1 1]);

r = wec1.Rcir;

wam_run = WamitRunCondition(path, name);           
wam_run.Rho = rho;                   

wam_run.FloatingBodies = wec1;               

T = IWaves.Lam2T(lam, h);
wam_run.T = T;
wam_run.Beta = Beta;                                   
wam_run.H = h;

r = r + dr*a;

z = -h*(1-cos(0:pi/2/nZ:pi/2));
z = round(z*10^6)/10^6;
points = makeCirWFPoints(r, nTheta, z);

wam_run.FieldPoints = points;   

wam_run.WriteRun;   

%% Run Wamit

wam_run.Run;

%% Read run f

wam_result = WamitResult(wam_run);                   
wam_result.ReadResult;                              

waveC = wam_result.WavePoints;
hydroF = wam_result.FreqDomForces;
floatB = wam_result.FloatingBodies;

save([mainPath 'cyl6dof_WamOutUT'], 'waveC', 'hydroF', 'floatB');

%% Create the cylinder hydrobody

load cyl6dof_WamOutUT;

%%

cylHydBod = computeHydroBody(waveC, hydroF, floatB, 'SigFigCutoff', 6, 'AccTrim');

mainPath = fileparts(which('mwave.home'));
mainPath = [mainPath '\UnitTests\files\'];

save([mainPath 'cyl6dof_hbUT'], 'cylHydBod');

%% 2 - Compute for single body verification

% Write WAMIT run
lam = 10;

Beta = linspace(0,pi/2,5);

mainPath = fileparts(which('mwave.home'));
mainPath = [mainPath '\UnitTests\files\'];

% Circ wave fields - write runs

% Cylinder
name = 'ver1';
path = [mainPath name];

h = 10;     
rho = 1000;

rad = 1;
hei = 1.5;
draft = 1;

Ntheta = 64;
Nr = 8;
Nz = 12;
    
wec1 = FloatingCylinder(rho, rad, hei, draft, Ntheta, Nr, Nz);

wec1.Handle = 'cyl';     
wec1.Modes = ModesOfMotion([1 1 1 1 1 1]);

r = wec1.Rcir;

wam_run = WamitRunCondition(path, name);           
wam_run.Rho = rho;                   

wam_run.FloatingBodies = wec1;               

T = IWaves.Lam2T(lam, h);
wam_run.T = T;
wam_run.Beta = Beta;                                   
wam_run.H = h;

wam_run.FieldArray = BemFieldArray([-20 -20 0], [0.4 0.4 1], [101 101 1]);

wam_run.WriteRun;   

%% Run Wamit

wam_run.Run;


%% Read run f

wam_result = WamitResult(wam_run);                   
wam_result.ReadResult;                              

wWave1cyl = wam_result.WaveArray;
hydroF = wam_result.FreqDomForces;
floatB = wam_result.FloatingBodies;

wComp1cyl = FreqDomComp(hydroF, floatB);

save([mainPath 'cyl1_verUT'], 'wComp1cyl', 'wWave1cyl');

%% 3 - Compute for two body verification

% Write WAMIT run
lam = 10;

Beta = linspace(0,pi/2,5);

mainPath = fileparts(which('mwave.home'));
mainPath = [mainPath '\UnitTests\files\'];

% Circ wave fields - write runs

% Cylinder
name = 'ver2';
path = [mainPath name];

h = 10;     
rho = 1000;

rad = 1;
hei = 1.5;
draft = 1;

Ntheta = 64;
Nr = 8;
Nz = 12;
    
wec1 = FloatingCylinder(rho, rad, hei, draft, Ntheta, Nr, Nz);
wec1.XYpos = [0 -4];
wec1.Handle = 'cyl';     
wec1.Modes = ModesOfMotion([1 1 1 1 1 1]);

wec2 = FloatingBody(wec1);
wec2.XYpos = [0 4];

wam_run = WamitRunCondition(path, name);           
wam_run.Rho = rho;                   

wam_run.FloatingBodies = [wec1 wec2];               

T = IWaves.Lam2T(lam, h);
wam_run.T = T;
wam_run.Beta = Beta;                                   
wam_run.H = h;

wam_run.FieldArray = BemFieldArray([-20 -20 0], [0.4 0.4 1], [101 101 1]);

wam_run.WriteRun;   

%% Run Wamit

wam_run.Run;


%% Read run f

wam_result = WamitResult(wam_run);                   
wam_result.ReadResult;                              

wWave2cyl = wam_result.WaveArray;
hydroF = wam_result.FreqDomForces;
floatB = wam_result.FloatingBodies;

wComp2cyl = FreqDomComp(hydroF, floatB);

save([mainPath 'cyl2_verUT'], 'wComp2cyl', 'wWave2cyl');




