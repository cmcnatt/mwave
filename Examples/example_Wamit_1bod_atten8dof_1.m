% This test sets up and runs Wamit on a body with 8 DoF using the Wamit
% generalized modes. It used a built-in mwave geometry which makes setting
% up the body relatively straightforward, and some amount of explanation is
% given on what goes on under the hood.
%
% For a more introductory case with more explanations, look at 
% example_Wamit_1bod_6dof_2

clear; close all; clc;

%% Set up the run

% Currently, this example is not working and should be fixed at some
% point... TODO: fix example
warning('Example known not to be working');

run_name = 'wam_1b_atten_1';         
folder = [mwavePath '\Examples\bemRunFolder'];  
if 0 == exist(folder, 'dir')
    mkdir(folder);
end
delete([folder '\*']);

rho = 1000;     

% Here we're creating a "Hinged Barge" which I generally also call an
% attenuator. This is the geometry (with different parameters) that I used
% in the Ocean Engineering paper 2014
len = 60;
dia = 10;
hingePos = [-10; 10];   % two hinges, these are the x positions. This 
                        % creates 3 bodies each with same length
sphereRad = 0.1*len;    % How large the spherical ends are

Nx = 80;                % number of panels along body length
Ntheta = 24;            % number of panels around the body

wec = FloatingSphereEndCyl(rho, len, dia/2, sphereRad, hingePos, ...
    Nx, Ntheta, 'notch');  
% What take a little while is computing the mass matrix. It is an 8x8 mass
% matrix and there is coupling between some of the modes. This is computed
% using the MassBody and MassPoint classes and following the method
% described in the Ocean Engineering paper. Also see the Newman reference
% from that paper. 

% Let's look at the Geometry..
figure;
plot(wec.PanelGeo);
axis equal;
title('Two-hinge Attenuator Panelization');
set(gca, 'view', [-35 15]);


% The mass matrix assumes the body is like a sausage, and what is show in
% the panelization is also mirrored above the surface for computing mass.

wec.Handle = 'atten';
wec.Modes = ModesOfMotion([1 1 1 1 1 1 1 1]);   % here, since the body has 
                                                % 8 DoF, we have to supply
                                                % a vector of length 8
                                                
% The extra modes of motion are handled by the Wamit Generalized modes
% functionality which relies on the Newmodes.dll - In this case,
% Newmodes.dll has been modified by me, so you need to use the same
% Newmodes.dll that is supplied with mwave. It needs to be copied into the
% same folder as the Wamit.exe

% The floating body has:
wec.WamIGenMds;
% which is just a number that is given to Wamit so that it knows which
% generalized modes function to use in Newmodes. 

% In this case, the FloatingSphereEndHinge body also creates another file
% called name_xhinge.dat, that gives the orientation of the body and the
% location of the hinges. This is read in by the newmodes function. It
% shouldn't be necessary except for a bug in Wamit.
                                            
wam_run = WamitRunCondition(folder, run_name);  

wam_run.Rho = rho;      
wam_run.T = 4:0.1:12;   
wam_run.Beta = 0;       
wam_run.H = Inf;        
wam_run.FloatingBodies = wec;       

wam_run.WriteRun;                   

%% Run Wamit

wam_run.Run;                           

%% Read results

wam_result = WamitResult(wam_run);  
wam_result.ReadResult;              

hydroForces = wam_result.FreqDomForces;
% Here, we have 8 DoF
A = hydroForces.A;      
B = hydroForces.B;      
C = hydroForces.C;      
Fex = hydroForces.Fex;  

T = hydroForces.T;     
figure;
subplot(3,1,1);
% DoF 7 and 8 are the hinge modes (or flex)
% Because of symmetry the radiation forces 7-7 and 8-8 are the same
axe = plotyy(T, [squeeze(A(:,7,7)) squeeze(A(:,8,8))], T, ...
    squeeze(A(:,3,7)));
title({'Results for Attenuator', 'Added Mass'});
ylabel(axe(1), 'kg*m^2');
ylabel(axe(2), 'kg*m');
legend('flex 1-flex 1', 'flex 2-flex 2', 'heave-flex 1');
subplot(3,1,2);
axe = plotyy(T, [squeeze(B(:,7,7)) squeeze(B(:,8,8))], T, ...
    squeeze(B(:,3,7)));
title('Hydrodynamic Damping');
ylabel(axe(1), 'Ns/m*m^2');
ylabel(axe(2), 'Ns/m*m');
legend('flex 1-flex 1', 'flex 2-flex 2', 'heave-flex 1');
subplot(3,1,3);
plot(T, abs([squeeze(Fex(:,1,7)) squeeze(Fex(:,1,8))]));
title('Excitation Force');
legend('flex 1', 'flex 2');
xlabel('Period (s)');
ylabel('Nm');


%% Analyze results


hydroComp = FreqDomComp(hydroForces, wec);

dof = wec.Modes.DoF;
Dpto = zeros(dof, dof);     
Dpto(7,7) = 10^8;           % In this case, we set modes 7 and 8 (the 
Dpto(8,8) = 10^8;           % hinges) to have a rotational PTO damping
           

hydroComp.SetDpto(Dpto);    
RAO = hydroComp.Motions;    
power = hydroComp.Power;    
                            
power1 = squeeze(power(:,1,7)); % power from just the first hinge
power2 = squeeze(power(:,1,8)); % power from the second hinge
powerT = power1 + power2;       % total power absorbed by the body

incFlux = IWaves.UnitEnergyFlux(rho, T, hydroComp.H).';    
RCW = powerT./incFlux./dia;

figure;
subplot(3,1,1);
axe = plotyy(T, abs(squeeze(RAO(:,1,3))), T, ...
    180/pi*abs([squeeze(RAO(:,1,7)) squeeze(RAO(:,1,8))]));
ylabel(axe(1), 'Heave (m/m)');
ylabel(axe(2), 'Flex (deg/m)');
title({'Results for Attenuator', 'Response Amplitude Operator'});
legend('Heave', 'Flex 1', 'Flex 2', 'NorthWest');

subplot(3,1,2);
plot(T, [power1, power2, powerT]./1000);
title('Power in 2 m high waves')
ylabel('kW');
legend('Hinge 1', 'Hinge 2', 'Total');

subplot(3,1,3);
plot(T, RCW);
title('Relative capture width (Total)');
xlabel('Period (s)');
ylabel('RCW');

