function [tank] = curvedTankPoints()

r = 9;          % radius of curved tank curve in meters
botLen = 12;    % length of the beach side of the curved tank
sidLen = 4.4;   % length of the glass wall side of the curved tank

wgcen = [-0.25 6];  % approximate (x,y) location of the center of the wave 
                    % gauge array relative the center of the tank curve

angStart = 64;  % approximate angle (in degrees) of the start of the curve 
                % of the wave tank

% set up tank with wave makers along top
thetaStart = pi/180*angStart;

topRight = r*[cos(thetaStart) sin(thetaStart)];
botRight = topRight - [0 sidLen];
botLeft = botRight - [botLen 0];

thetaStop = atan2(botLeft(2), botLeft(1));
theta = linspace(thetaStart, thetaStop, 101)';

curve = r*[cos(theta) sin(theta)];

tank = [curve; botRight; curve(1,:)];
np = size(tank,1);

% center about wave gauges
tank = tank - ones(np,1)*wgcen;

% rotate the tank
R = [0 -1; 1 0];    % 90 deg rotation matrix;

for n = 1:np
    tank(n,:) = (R*tank(n,:)')';
end

end