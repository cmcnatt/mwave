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
% This test sets up and runs Nemoh on a flap and plot the wave field. 
%
% For more information, see
%   - Wamit_1bod_6dof_1
%   - Wamit_1bod_6dof_2
%   - Nemoh_1bod_6dof_1
%   - Nemoh_createHB_1

%% Set up the run
        
folder = [mwavePath 'Examples\BemRuns\'];  
                                                    
rho = 1000;     
h = 50;

len = 5;
beam = 30;
draft = 20;

Nx = 5;            
Ny = 30;                 
Nz = 20;                
flap = FloatingFlap(rho, len, beam, draft, Nx, Ny, Nz);
flap.Modes = ModesOfMotion([0 0 0 0 1 0]);

nps = [101 101 1];
lens = [100 100 1];
dels = lens./(nps - [1 1 0]);
starts = -lens./2;
starts(3) = 0;

fieldArray = BemFieldArray(starts, dels, nps);


%% Run Nemoh

nem_name = 'nem_1b_1dof_1'; 

% Create a NemohRunCondition object - this builds the Wamit run
nem_run = NemohRunCondition([folder nem_name]); % The NemohRunCondition 
                                                % just takes a folder that 
                                                % is the exact location 
                                                % where the run goes
                                                
nem_run.ComputeHS = false;
nem_run.Rho = rho;      
nem_run.H = h;       

nem_run.OmegaLims = [0.5, 1.5]; % about 12-4 s wave periods
nem_run.OmegaCount = 3;        

% The wave directions are the same as the frequency
nem_run.BetaLims = [0 pi/4];       
nem_run.BetaCount = 2;

nem_run.FloatingBodies = flap;       

nem_run.FieldArray = fieldArray;

% Set the location of the Nemoh executables
nem_run.ExePath = 'C:\nemoh';

nem_run.WriteRun;           

nem_run.Run       

nem_result = NemohResult(nem_run);  
nem_result.ReadResult;          

hfn = nem_result.FreqDomForces;
wfn = nem_result.WaveArray;

%% Run Wamit - for comparison

wam_name = 'nem_1b_1dof_1w'; 

wam_run = WamitRunCondition([folder wam_name], wam_name); 
wam_run.Rho = nem_run.Rho;      
wam_run.H = nem_run.H;       

wam_run.T = nem_run.T;
wam_run.Beta = nem_run.Beta;

wam_run.FloatingBodies = nem_run.FloatingBodies;

wam_run.FieldArray = nem_run.FieldArray;

wam_run.WriteRun;           

wam_run.Run       

wam_result = WamitResult(wam_run);  
wam_result.ReadResult;          

hfw = wam_result.FreqDomForces;
wfw = wam_result.WaveArray;

%% Compare results

[X, Y] = wfn.FieldPoints;
omega = 2*pi./wfn.T;
beta = hfn.Beta;

etan = wfn.Elevation('Diffracted');
etaw = wfw.Elevation('Diffracted');

figure;
for m = 1:3
    for n = 1:2
        rowind = 2*(2*(m-1)+n-1);
        subplot(6,2,rowind + 1)
        pcolor(X, Y, real(etan{m, n}));
        fet;
        set(gca, 'xtick', [], 'ytick', []);
        ylabel({['\omega = ' num2str(omega(m), '%4.2f')], ['\beta = ' num2str(beta(n), '%4.2f')]})
        if (rowind == 0)
            title('Nemoh');
        elseif (rowind == 10)
            set(gca, 'xtick', [-50 0 50], 'ytick', [-50 0 50]);
        end
        subplot(6,2,rowind + 2)
        pcolor(X, Y, real(etaw{m, n}));
        fet;
        set(gca, 'xtick', [], 'ytick', []);
        if (rowind == 0)
            title('Wamit');
        end
    end
end

etan = wfn.Elevation('Radiated');
etaw = wfw.Elevation('Radiated');

figure;
for m = 1:3
    rowind = 2*(m-1);
    subplot(3,2,rowind + 1)
    pcolor(X, Y, real(etan{m}));
    fet;
    set(gca, 'xtick', [], 'ytick', []);
    ylabel({['\omega = ' num2str(omega(m), '%4.2f')]})
    if (rowind == 0)
        title('Nemoh');
    elseif (rowind == 10)
        set(gca, 'xtick', [-50 0 50], 'ytick', [-50 0 50]);
    end
    subplot(3,2,rowind + 2)
    pcolor(X, Y, real(etaw{m}));
    fet;
    set(gca, 'xtick', [], 'ytick', []);
    if (rowind == 0)
        title('Wamit');
    end
end