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
% In this example, we'll use the HydroBody we created in example
% Wamit_createHB_1 to create an array and compute the relavent hydrodynamic
% forces, power, and the wave field. 
%
% For more about a HydroBody, see hydroBody_1

%% Load the HydroBody

% load 'wam_hb1_1_hb';    % Load the HydroBody based on its name..

run_name = 'wam_hb_1';  % Or we can be more explicit if necessary         
folder = [mwavePath '\Examples\WamitRuns\' run_name];  
load([folder '\wam_hb1_1_hb']);

%% Create an array

% First we need to figure out the spatial arrangement of the array. Let's
% do 5 bodies arranged like this
%
%   *(-50,100)
%
%               *(50,50)
%
%   *(-50,0)
%
%               *(50,-50)
%
%   *(-50,-100)

% We create an array of HydroBodies 
for m = 1:5
    hbs(m) = HydroBody(hydBody);    % Here we use the copy constructor to 
                                    % make new HydroBodies that are copies
                                    % of the HydroBody that we just loaded.
                                    % If we said: hbs(n) = hydBody;
                                    % everytime we changes one of the
                                    % bodies, it would change them all
end

% Then we set their positions
hbs(1).XYpos = [-50, 100];
hbs(2).XYpos = [-50, 0];
hbs(3).XYpos = [-50, -100];
hbs(4).XYpos = [50, 50];
hbs(5).XYpos = [50, -50];

% And that's our array of bodies..

%% Create array computation object

% To make a computation, the array needs some incident waves. For more
% info, see the example hydroBody_1
T = hydBody.T;      
h = hydBody.H;      
beta = [0 pi/4 pi/3, pi/2]; 
a = ones(size(T)); % Unit amplitude waves
for m = 1:length(beta)
    iwaves(m) = PlaneWaves(a, T, beta(m), h);
end

% HydroArrayComp, like HydroBodyComp does all of the heavy lifting. 
arrayComp = HydroArrayComp(hbs, iwaves);    % It takes an array of 
                                            % HydroBodies and the incident 
                                            % waves in its constructor

% HydroArrayComp does not actually compute any values until you ask it for
% something. So, let's do that..

%% Compute array values

% This will trigger the computation. Unless the incident waves or the
% bodies are changed, this does not have to be redone:
A = arrayComp.A;        % The added mass matrix for the entire array, 
                        % (nT x DoF x DoF), where DoF in this case is 
                        % 40 = 5 bodies * 8 DoF each
B = arrayComp.B;        % Hydrodynamic damping for whole array
C = arrayComp.C;        % Hydrostatic stiffness 
Fex = arrayComp.Fex;    % Excitation force for the entire array, 
                        % (nT x nI x DoF), nI is the number of incident
                        % waves

% You may get a warning message:
%
% Warning: Matrix is close to singular or badly scaled. Results may be 
% inaccurate. RCOND =  6.653070e-65. 
%
% This is because the diffraction transfer matrix has some small values,
% which are likely meaningless, and which makes the matrix solution hard to
% find. To get rid of these values, when you create the HydroBody, use, the
% options 'SigFigCutoff', 5, and 'AccTrim'. See example, Wamit_createHB_1

% Compute power
% set the PTO damping
d = 10^8;           % damping value
dof = hydBody.DoF;
Dpto = zeros(dof, dof);     
for m = 0:4         % Loop through to set the damping on each body
    ih1 = m*8+7;    % index of hinge 1 for a given body
    ih2 = m*8+8;    % index of hinge 2
    Dpto(ih1,ih1) = d;           
    Dpto(ih2,ih2) = d;           
end
arrayComp.SetDpto(Dpto);    

power = arrayComp.Power;    
ihinge1 = 7:8:dof;  % another way to index the hinge modes
ihinge2 = 8:8:dof;
% Get the total power in the array for each wave period, direction
power = squeeze(power(:,:,ihinge1)) + squeeze(power(:,:,ihinge2));   

figure;
plot(T, [squeeze(power(:,1)), squeeze(power(:,2)), squeeze(power(:,3)), ...
    squeeze(power(:,4))]./1000);
title('Power in 2 m high waves')
ylabel('kW');
xlabel('Period (s)');
Betas = {'\beta = 0', '\beta = \pi/4', '\beta = \pi/3', '\beta = \pi/2'};
legend(Betas);

%% Create wave fields. 
isarray = true;
x = -200:2:200;
[X, Y] = meshgrid(x, x);    % square

waveField = arrayComp.WaveField(isarray, X, Y, 'NoVel');    

%%
eta = waveField.Elevation('Total');             

sect = hydBody.WaterPlaneSec;
onesSect = ones(size(sect, 1), 1);

thet = (0:(2*pi/50):2*pi)';
cir = hydBody.Rcir*[cos(thet) sin(thet)];
onesCir = ones(size(cir, 1), 1);

iT = 4;
figure;
for m = 1:4
    subplot(2,2,m);
    pcolor(X,Y,abs(eta{iT,m}));
    hold on;
    for n = 1:5
        % The sections are in body coordinates and so need to me moved to
        % where the body is in global coordinates
        thisSect = sect + onesSect*hbs(n).XYpos;    
        thisCir = cir + onesCir*hbs(n).XYpos;
        plot(thisSect(:,1), thisSect(:,2), 'w');
        plot(thisCir(:,1), thisCir(:,2));
    end
    fet;
    set(gca, 'clim', [0.7 1.3]);
    title(Betas{m});
end

