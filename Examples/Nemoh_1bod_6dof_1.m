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
% This test sets up and runs Nemoh on a 6 dof cylinder. 
% The cylinder's geometry file is made in matlab with mwave. 
%
% For more information, see
%   - Wamit_1bod_6dof_1
%   - Wamit_1bod_6dof_2
%   - Nemoh_createHB_1

%% Set up the run
        
folder = [mwavePath 'Examples\BemRuns\'];  
                                                    
rho = 1000;    
h = 50;

diameter = 10;
draft = 10;
height = 15;            

Ntheta = 48;            
Nr = 8;                 
Nz = 24;                

cyl = FloatingCylinder(rho, diameter/2, height, draft, Ntheta, Nr, Nz); 

cyl.Handle = 'cyl';

%% Run Nemoh

nem_name = 'nem_1b_6dof_1'; 

% Create a NemohRunCondition object - this builds the Wamit run
nem_run = NemohRunCondition([folder nem_name]); % The NemohRunCondition 
                                                % just takes a folder that 
                                                % is the exact location 
                                                % where the run goes
                                                
nem_run.ComputeHS = false;
nem_run.Rho = rho;      
nem_run.H = h;       

% The NemohRunCondition creates a Nemoh.cal file, which is used by the
% preProcessor. The Nemoh.cal file takes radial frequency limits and count
nem_run.OmegaLims = [0.5, 1.5]; % about 12-4 s wave periods
nem_run.OmegaCount = 24;        % 24 frequencies - evenly spaced in omega

% The wave directions are the same as the frequency
nem_run.BetaLims = [0 0];       % single direction of 0 degrees
nem_run.BetaCount = 1;

nem_run.FloatingBodies = cyl;       

% Set the location of the Nemoh executables
nem_run.ExePath = 'C:\nemoh';

nem_run.WriteRun;           

%%
tic
nem_run.Run
nemRT = toc;

nem_result = NemohResult(nem_run);  
nem_result.ReadResult;          

hfn = nem_result.FreqDomForces;

%% Run Wamit - for comparison

wam_name = 'nem_1b_6dof_1w'; 

wam_run = WamitRunCondition([folder wam_name], wam_name); 
wam_run.Rho = nem_run.Rho;      
wam_run.H = nem_run.H;       

wam_run.T = nem_run.T;
wam_run.Beta = nem_run.Beta;

wam_run.FloatingBodies = nem_run.FloatingBodies;

wam_run.WriteRun;           

wam_run.Run       

wam_result = WamitResult(wam_run);  
wam_result.ReadResult;          

hfw = wam_result.FreqDomForces;

%% Compare results

omega = 2*pi./hfn.T;

An = hfn.A;
Bn = hfn.B;
Cn = hfn.C;
Fexn = hfn.Fex;

Aw = hfw.A;
Bw = hfw.B;
Cw = hfw.C;
Fexw = hfw.Fex;

figure;

inds = [1 1; 2 2; 3 3; 4 4; 5 5; 6 6];
titles = {'Surge-surge', 'Sway-sway', 'Heave-heave', 'Roll-roll', ...
    'Pitch-pitch', 'Yaw-yaw'};

for n = 1:6
    
    i1 = inds(n, 1);
    i2 = inds(n, 2);

    subplot(6,2,2*(n-1)+1);
    plot(omega, [squeeze(An(:,i1,i2)) squeeze(Aw(:,i1,i2))]);
    ylabel(titles{n});
    set(gca, 'xtick', []);
    if (n == 1)
        title('Added Mass');
        legend('Nemoh', 'Wamit');
    elseif (n == 6)
        xlabel('\omega (rad/s)');  
        set(gca, 'xtick', (0.5:0.2:1.5));
    end
    subplot(6,2,2*(n-1)+2);
    plot(omega, [squeeze(Bn(:,i1,i2)) squeeze(Bw(:,i1,i2))]);
    set(gca, 'xtick', []);
    if (n == 1)
        title('Damping');
    elseif (n == 6)
        xlabel('\omega (rad/s)');   
        set(gca, 'xtick', (0.5:0.2:1.5));
    end
end

figure;

inds = [1 1; 2 2; 3 3; 4 4; 5 5; 6 6];
titles = {'Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'};

for n = 1:6
    
    subplot(6,1,n);
    plot(omega, abs([squeeze(Fexn(:,1,n)) squeeze(Fexw(:,1,n))]));
    ylabel(titles{n});
    set(gca, 'xtick', []);
    if (n == 1)
        title('Excitation Force');
        legend('Nemoh', 'Wamit');
    elseif (n == 6)
        xlabel('\omega (rad/s)');  
        set(gca, 'xtick', (0.5:0.2:1.5));
    end
end
