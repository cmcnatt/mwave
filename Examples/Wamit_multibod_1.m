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
% This example sets up and computes an array in Wamit directly - that is,
% we put multiple bodies into Wamit and Wamit computes the full Added Mass,
% Damping, and Excitation force matrices. This is as oppposed to using
% Interation Theory.
%
% For a more introductory case with more explanations, look at 
% Wamit_1bod_6dof_2
%
% For more introduction on WaveFields take a look at
% Wamit_1bod_wavefield_1

%% Set up run

run_name = 'wam_mb_1';         
folder = [mwavePath '\Examples\BemRuns\' run_name];  

rho = 1000;     

% We'll have 3 bodies in the run - 2 attenuators and 1 flap. The
% attenuators will be the same, except that one will be turned
len = 60;
dia = 10;
hingePos = [-10; 10];   
sphereRad = 0.1*len;    

Nx = 80;                
Ntheta = 24;            

atten1 = FloatingSphereEndCylHinge(rho, len, dia/2, sphereRad, hingePos,...
    Nx, Ntheta);  
atten1.Modes = ModesOfMotion([1 1 1 1 1 1 1 1]);   
atten1.Handle = 'atten1';

% Here we want to make another attenuator exactly like the first, except
% rather than using 3 lines, we'll use the copy constructor in FloatingBody
atten2 = FloatingBody(atten1);
atten2.Handle = 'atten2';

% If you just said:
% atten2 = atten1;
% that would create a shallow copy and changes made to atten2 would change
% atten1

% Place the bodies in (x,y) global coordinates. Currently they are both at
% (0,0), which isn't going to work

atten1.XYpos = [-40 40];    % This is whay .XYpos is separated from .Zpos
atten2.XYpos = [0 -40];

% And let's rotation atten1 by -45 degrees (i.e. rotate it clockwise)
atten1.Angle = 315;          % degrees!  This probably needs to be changed

% And we'll throw a cylinder into the mix
diameter = 10;
draft = 10;
height = 15;            

Ntheta = 48;            
Nr = 8;                 
Nz = 24;                

cyl = FloatingCylinder(rho, diameter/2, height, draft, Ntheta, Nr, Nz);
cyl.Handle = 'cylinder';
cyl.Modes = ModesOfMotion([1 1 1 1 1 1]);
cyl.XYpos = [10 0];

% Now we create an array of the floating bodies to pass to
% WamitRunCondition
bodies = [atten1, atten2, cyl];

wam_run = WamitRunCondition(folder, run_name);  

wam_run.Rho = rho;      
wam_run.T = 7;      % Just do a single period (same some time)   
wam_run.Beta = 0;       
wam_run.H = 100;        

% Here we tell Wamit about our array.
wam_run.FloatingBodies = bodies;     

% Let's also create a wave field - makes it nice to see the bodies
wam_run.FieldArray = ...
    BemFieldArray([-100 -100 0], [1 1 1], [201 201 1]);

wam_run.WriteRun;                   

%% Run Wamit

%wam_run.Run;                           
 wam_run.Run('Background');   
 
%% Read results

wam_result = WamitResult(wam_run);  
wam_result.ReadResult;              

hydroForces = wam_result.HydroForces;   % The hydroForces has 22 DoF - 8 
                                        % for each attenuator and 6 for the
                                        % cylinder
waveField = wam_result.WaveArray;   

%% Analyze Results

% In this case, we have to give the HydroBodyComp class the array of
% FloatingBodies
hydroComp = HydroBodyComp(hydroForces, bodies);

dof = hydroComp.DoF;
dpto = 10^8;
Dpto = zeros(dof, dof);     % Again, there are 22 DoF     
Dpto(7,7) = dpto;           % Set PTO on atten1 
Dpto(8,8) = dpto;
Dpto(15,15) = dpto;         % Mode 15 is the first hinge of the second 
                            % attenuator. 8 dof for first atten, then the
                            % hinge is the 7th mode of the second
Dpto(16,16) = dpto;
% We'll leave the cylinder's pto damping alone
hydroComp.SetDpto(Dpto);

% Since we just have one period - we'll write the power output to the
% screen
power = hydroComp.Power;    
fprintf('\nPower from attenuators:\n')
fprintf('Atten 1, front hinge: %3.0f kW\n', power(1,1,7)./1000);
fprintf('Atten 1, back hinge: %3.0f kW\n', power(1,1,8)./1000);
fprintf('Atten 2, front hinge: %3.0f kW\n', power(1,1,15)./1000);
fprintf('Atten 2, back hinge: %3.0f kW\n\n', power(1,1,16)./1000);


[X, Y] = waveField.FieldPoints;  
waveField.BodyMotions = hydroComp.Motions;

etaT = waveField.Elevation('Total');

figure;
pcolor(X,Y,abs(etaT{1}));
fet;
set(gca,'clim',[0.6 1.4]);
title({'Wamit array computation', 'Abs Wave elevation'});
xlabel('x (m)');
ylabel('y (m)');
                                    
