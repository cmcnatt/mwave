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
% In this run, we will create a 'HydroBody' in Wamit. A HydroBody is a
% special class that is required for the array interaction theory. It is
% basically a FloatingBody with hydrodynamic parameters attached, including
% added mass, damping, hydrostatics, and the so-called "force transfer
% matrix" which is described in the 2014 Ocean Engineering paper. In
% addition to these it has wave field coefficient - Radiation coefficients
% and the diffraction transfer matrix. HydroBodies depend on wave period
% and water depth (and geometry), but not incident wave direction
%
% It is created by computing the wavefield on a circular cylinder of points
% circumscribing a geometry, and then running waves at the geometry from
% all directions surrounding the body. At this point, symmetry of the body
% is not exploited, but it should be. This is a big TODO
%
% For more info on the theory see the 2013 IJOME paper and the 2014 Ocean
% Engineering paper.
%
% For a more introductory case with more explanations, look at 
% Wamit_1bod_6dof_2
%
% For more introduction on WaveFields take a look at
% Wamit_1bod_wavefield_1
%
% For more info on the geometry used here, check out Wamit_1bod_atten8dof_1

%% Set up run

run_name = 'wam_hb_1';         
folder = [mwavePath '\Examples\BemRuns\' run_name];  

rho = 1000;     

% set up a normal body
len = 60;
dia = 10;
hingePos = [-10; 10];   
sphereRad = 0.1*len;    

Nx = 80;                
Ntheta = 24;            

wec = FloatingSphereEndCylHinge(rho, len, dia/2, sphereRad, hingePos, ...
    Nx, Ntheta);  

% Now we need to give Wamit lots of incident wave directions covering the
% full range from 0 to 2pi (Like I said in the intro, body geometry
% symmetry has not been taken advantage of yet. That's a TODO)
N = 40;                         % 40 directions
Beta = 0:2*pi/N:2*pi*(1-1/N);   % from 0 to 2pi

% Now we have to get the wave field at points circumscribing the body
r = wec.Rcir;                   % Rcic is the radius of the circle that 
                                % circumscribes the body
dr = 0.1;                       % And we need the radius to be slightly 
r = r + dr;                     % bigger than that, because you don't want
                                % to compute field points too close to the
                                % body in Wamit                   
nZ = 200;                       % The number of points in the z direction
nTheta = 2^8;                   % The number of points in the angular dir

h = 50;                         % water depth (m). We need to set it here 
                                % because we need it to create field points
                                % CANNOT DO INFITE DEPTH - the theory is
                                % not set up to support the Inf depth case
                                % yet, but that is a TODO

% Now we actually creat the (x,y,z) points of the cylindrical surface
z = -h*(1-cos(0:pi/2/nZ:pi/2)); % Use cosine spacing in z, because you want
                                % more points near the free-surface because
                                % progressive waves decay exponentially
                                % with depth - see 2013 IJOME paper
z = round(z*10^6)/10^6;         % round so that the top point is at z = 0
points = makeCirWFPoints(r, nTheta, z); % This function makes the points. 
                                        % Only hard because they have to be
                                        % in one long Nx3 list

wam_run = WamitRunCondition(folder, run_name);  

wam_run.Rho = rho;      
wam_run.T = 4:1:12;             % Do fewer periods   
wam_run.Beta = Beta;            % The list of incident directions       
wam_run.H = h;        
wam_run.FloatingBodies = wec;       

wam_run.FieldPoints = points;   % Here is where we actually apply the field
                                % points we created to the Wamit run
wam_run.WriteRun;               % Good to go

%% Run Wamit

% wam_run.Run;                           
wam_run.Run('Background');           

%% Read results, and save useful objects

wam_result = WamitResult(wam_run);  
wam_result.ReadResult;          % Reading the results can take a while 
                                % because they are massive text files, and
                                % they take up a lot space: ~GB

waveCir = wam_result.WavePoints;
hydroForces = wam_result.HydroForces;

% It can be a good idea to just say these, so that they can be reloaded 
% later, but they are also massive files.. 
% save([folder '\wam_hb1_1_out'], 'waveCir', 'hydroForces', 'wec');

%% create HydroBody and save it

% load 'wam_hb1_1_out'

% This is the function that does it all.. It has a number of different
% options and really needs to be documented a bit better..
hydBody = computeHydroBody(waveCir, hydroForces, wec, 'SigFigCutoff', 5,...
    'AccTrim');
% 'SigFigCutoff', N tells the computation to only use the first N 
% significant figures. 'AccTrim' removes values on the ends that do not
% contribute significantly to reproducing the wave field. Without these
% options, when you compute arrays, you may get a warning message:
%
% Warning: Matrix is close to singular or badly scaled. Results may be 
% inaccurate. RCOND =  6.653070e-65. 
%
% This is because the diffraction transfer matrix has some small values,
% which are likely meaningless, and which makes the matrix solution hard to
% find. To get rid of these values, when you create the HydroBody, use, the
% options 'SigFigCutoff', 5, and 'AccTrim'. See example, Wamit_createHB_1


% Here, we're saving it in the working folder that we've been using and
% saving it by a name that will help us recognize it. It does take a while
% to do these runs in Wamit, so you really only want to create a HydroBody
% once and then save it 
save([folder '\wam_hb1_1_hb'], 'hydBody');

% And that's it for this example, in other examples we'll use the HydroBody
% to see its power!!
%   - hydroBody_1: use a single HydroBody to compute power and WaveFields
%   - array_inter_1: use the HydroBody to compute array interactions