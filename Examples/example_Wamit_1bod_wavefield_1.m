% This test sets up and runs Wamit on a 6 dof cylinder. 
% In this case, we want to compute the wave field around the cylinder. 
% For a more introductory case with more explanations, look at 
% example_Wamit_1bod_6dof_2

clear; close all; clc;

%% Set up the run

run_name = 'wam_1b_wf_1';         
folder = [mwavePath '\Examples\bemRunFolder'];  
if 0 == exist(folder, 'dir')
    mkdir(folder);
end
delete([folder '\*']);

rho = 1000;     

% The cylinder parameters are:
diameter = 10;
draft = 10;
height = 15;            

Ntheta = 48;            
Nr = 8;                 
Nz = 24;                

cyl = FloatingCylinder(rho, diameter/2, height, draft, Ntheta, Nr, Nz);
cyl.Handle = 'Cylinder';                                      
cyl.Modes = ModesOfMotion([1 0 1 0 1 0]);   % No need to compute the modes
                                            % orthogonal to the indicident
                                            % wave direction
                                            
wam_run = WamitRunCondition(folder, run_name);  

wam_run.Rho = rho;      
wam_run.T = 4:1:12;        
wam_run.Beta = 0;       
wam_run.H = Inf;        
                        
wam_run.FloatingBodies = cyl;       

% There are two types of wave fields that can be measured with Wamit:
% 1) an array of points - this could be a three dimensional array, but
% mwave is only set up to handle 2D (x,y) arrays currently.
% 2) A list of (x,y,z) points. This exmple just looks at the array.
start = [-40 -40 0];    % start location, [x, y, z]
delta = [0.8 0.8 1];    % delta (point spacing) in [x, y, z]
nPoint = [101 101 1];   % Number of points in [x, y, z]. Having only 1  
                        % point in the z direction creates a 2D array
wam_run.FieldArray = BemFieldArray(start, delta, nPoint);
wam_run.ComputeVelocity = true; % Here we're telling Wamit to compute the 
                                % velocity at the field points. The default
                                % is not to. It does make it take longer,
                                % particularly with reading and writig the
                                % values.

wam_run.WriteRun;                   
                                    
%% Run Wamit

wam_run.Run;
% wam_run.Run('Background');             

%% Read results

wam_result = WamitResult(wam_run);  
wam_result.ReadResult;              % Reading wave field results can take a
                                    % a long time, because they may contain
                                    % lots of points.

hydroForces = wam_result.FreqDomForces;

% This is the wave field that results from the Wamit computation. 
waveField = wam_result.WaveArray;   

% It is a 'FBWaveField' for FloatingBody wave field. There are other more
% basic types of wavefields, including the "interface" IWaveField (really
% an abstract class), and the WaveField class. FBWaveFields distinguish
% themselves by have difference "sub" WaveFields, including the diffracted
% WaveField and radiated WaveFields in each DoF. These are the wave fields
% that are computed by Wamit. Wamit only computes, pressure and velocity. 
% Other wave field values are derived. 

% For more on basic WaveFields see the example waveFields_1

%% Analyze results (wave fields, no motions)

[X, Y] = waveField.FieldPoints;     % These are the locations where the 
                                    % wave field values were calculated. 
                                    % For an array it is a meshgrid.
                                    
waveField.BodyMotions = [];         % This sets the wave field's body
                                    % motions to be empty. It's not really
                                    % necessary because when you read the
                                    % wave field out of the WamitResutls
                                    % they are empty. But if you were to
                                    % run this section twice, the second
                                    % time around, the motions would have
                                    % been set and results as laid out
                                    % here wouldn't make sense.
                            
% The wave field has different types of 'Wave fields'. The Elevation method
% returns complex wave elevation values at the meshgrid of spatial
% locations. The outputs are cell arrays, where in a given cell is a given
% result (meshgrid of complex elevation values).
etaI = waveField.Elevation('Incident');     % incident (Nper x Ndir)
etaD = waveField.Elevation('Diffracted');   % diffracted (Nper x Ndir)
etaS = waveField.Elevation('Scattered');    % scattered, (Nper x Ndir) 
                                            % which is diff minus inc
etaR = waveField.Elevation('Radiated');     % radiated (Nper x DoF)

% The wave field can also compute the 'Total' wave field, but to do so, it
% needs to know what the body motions are. This will be shown below. 

% Since the motions are not set, the radiated wave fields are for
% unit-amplitude, zero-phase motion. 

% Let's look at the wave fields. We'll just look at 3 of the periods, and 2
% of the radiated wave fields (surge and heave)

iper = [1 4 9];
T = waveField.T;

showAbs = true;     % Look at the absolute value (abs), or look at the 
                    % real part (real) by setting showAbs to false

% Plot em
figure;
for m = 1:3
    subplot(5,3,m);
    if (showAbs)
        pcolor(X,Y,abs(etaI{iper(m)}));     % the incident waves look a 
                                            % little weird plotted like 
                                            % this because the amplitude is
                                            % 1 everwhere.
    else
        pcolor(X,Y,real(etaI{iper(m)}));
    end
    axis equal;
    axis tight;
    shading flat;
    colorbar;
    if (m == 2)
        title({'Wave fields without body motions',['Period = ' num2str(T(iper(m))) ' s']});
    else
        title(['Period = ' num2str(T(iper(m))) ' s']);
    end
    if (m == 1)
        ylabel('Incident');
    end
    
    subplot(5,3,3+m);
    if (showAbs)
        pcolor(X,Y,abs(etaD{iper(m)}));
    else
        pcolor(X,Y,real(etaD{iper(m)}));
    end
    axis equal;
    axis tight;
    shading flat;
    colorbar;
    if (m == 1)
        ylabel('Diffracted');
    end
    
    subplot(5,3,6+m);
    if (showAbs)
        pcolor(X,Y,abs(etaS{iper(m)}));
    else
        pcolor(X,Y,real(etaS{iper(m)}));
    end
    axis equal;
    axis tight;
    shading flat;
    colorbar;
    if (m == 1)
        ylabel('Scattered');
    end
    
    subplot(5,3,9+m);
    if (showAbs)
        pcolor(X,Y,abs(etaR{iper(m),1}));
    else
        pcolor(X,Y,real(etaR{iper(m),1}));
    end
    axis equal;
    axis tight;
    shading flat;
    colorbar;
    if (m == 1)
        ylabel('Radiated - Surge');
    end
    
    subplot(5,3,12+m);
    if (showAbs)
        pcolor(X,Y,abs(etaR{iper(m),2}));
    else
        pcolor(X,Y,real(etaR{iper(m),2}));
    end
    axis equal;
    axis tight;
    shading flat;
    colorbar;
    if (m == 1)
        ylabel({'Radiated - Heave', 'y (m)'});
        xlabel('x (m)')
    end
end

%% Analyze results (wave fields with motions)

% Now let's apply some real motions to this body. To compute the motions,
% use the FreqDomComp
hydroComp = FreqDomComp(hydroForces, cyl);

dof = cyl.Modes.DoF;
Dpto = zeros(dof, dof);     
Dpto(2,2) = 10^3;           % Set a damping in the heave mode of motion, 
                            % which is the second degree of freedom in this
                            % case

hydroComp.SetDpto(Dpto);   

% Set the wave field body's motion to the motions computed by hydroComp
xi = hydroComp.Motions;
waveField.BodyMotions = xi;

% Look at surge and heave
figure;
plot(T, abs(squeeze(xi(:,1,1:2))));
title('Body Motions');
xlabel('T (s)');
ylabel('RAO');
legend('Surge', 'Heave');

% Now, let's recomputed the radiated and compute the total. The incident, 
% diffracted, and scattered wave fields do not depend on body motions and 
% so will be the same, (and won't be plotted)
etaR = waveField.Elevation('Radiated');     % radiated (Nper x Ndir x DoF)
                                            % The radiated wave field
                                            % includes the waves due to
                                            % each actual body motion,
                                            % which includes direction
etaT = waveField.Elevation('Total');        % total (Nper x Ndir). Now we 
                                            % can also compute the total
                                            % wave field 
% Now let's look at them. 
figure;
for m = 1:3
    subplot(3,3,m);
    if (showAbs)
        pcolor(X,Y,abs(etaR{iper(m),1,1}));
    else
        pcolor(X,Y,real(etaR{iper(m),1,1}));
    end
    axis equal;
    axis tight;
    shading flat;
    colorbar;
    if (m == 2)
        title({'Wave fields with body motions',['Period = ' num2str(T(iper(m))) ' s']});
    else
        title(['Period = ' num2str(T(iper(m))) ' s']);
    end
    if (m == 1)
        ylabel('Radiated - Surge');
    end    

    
    subplot(3,3,3+m);
    if (showAbs)
        pcolor(X,Y,abs(etaR{iper(m),1,2}));
    else
        pcolor(X,Y,real(etaR{iper(m),1,2}));
    end
    axis equal;
    axis tight;
    shading flat;
    colorbar;
    if (m == 1)
        ylabel('Radiated - Heave');
    end
    
    subplot(3,3,6+m);
    if (showAbs)
        pcolor(X,Y,abs(etaT{iper(m)}));
    else
        pcolor(X,Y,real(etaT{iper(m)}));
    end
    axis equal;
    axis tight;
    shading flat;
    colorbar;
    
    if (m == 1)
        ylabel({'Total', 'y (m)'});
        xlabel('x (m)')
    end
end

%% Analyze Results (velocity)

% The wave field can also return pressure and velocity. For example
pT = waveField.Pressure('Total');       % Pressure in this case, is just 
                                        % rho*g*elevation
velT = waveField.Velocity('Total');     % Velocity is of size 
                                        % (Nper x Ndir x 3), where the 3 is
                                        % represents the x,y,z components.
                                        
inds = 1:5:101;        % Don't want to look at all the points

iper = 1;
U = velT{iper,1}(inds, inds);
V = velT{iper,2}(inds, inds);

figure;
quiver(X(inds, inds),Y(inds, inds),real(U), real(V));
axis equal;
axis tight;
title({'Velocity of Total Wave Field', ...
    ['Period = ' num2str(T(iper)) ' s']});
xlabel('x (m)');
ylabel('y (m)');

% It is also useful here to show the body in the wave field plot
sect = cyl.WaterPlaneSec;       % The water plane section is the x,y 
                                % coordinates of the body section that
                                % intersects with z = 0;
hold on;
plot(sect(:,1), sect(:,2),'k'); % You can see inside the section, the 
                                % velocity is zero

%% TODO: Analyze results (spectrum)

% You can also use the wave field as an irregular wave field

