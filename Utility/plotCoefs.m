function [] = plotCoefs(M, a, showAxis, varargin)

Na = length(a);
if ((Na-1)/2 ~= M)
    error('Wrong number of coeffs');
end

A = abs(a);
phi = angle(a);

for m = -M:M
    im = m + M + 1;
        
    x = [0 A(im)*sin(phi(im))];
    y = [m m];
    z = [0 A(im)*cos(phi(im))];
    
    plot3(x,y,z,varargin{:});
    hold on;
    plot3(x(2),y(2),z(2),'Marker', '.','MarkerSize', 10,varargin{:});
end

axis equal
axis off

 maxA = max(abs(A));

view([-45, 20]) 

if (showAxis)
    wid = 0.005;

    mArrow3([0 0 0], [0, -M-1, 0], 'stemWidth', wid);
    mArrow3([0 0 0], [0, M+1, 0], 'stemWidth', wid);
    mArrow3([0 0 0], [0, 0, 1.5*maxA], 'stemWidth', wid);

    thetaS = pi/2-pi/2/6;
    theta = linspace(0,thetaS,101);
    thetaS = thetaS+pi/2;
    theta = theta + pi/2*ones(size(theta));
    rad = 0.5*maxA;
    ys = zeros(size(theta));
    xs = rad*cos(theta);
    zs = rad*sin(theta);
    plot3(xs,ys,zs,'k','LineWidth', 50*wid);
    arlen = 0.25;
    ctS = cos(thetaS);
    stS = sin(thetaS);
    mArrow3(rad*[ctS 0 stS],  rad*[ctS-arlen*stS, 0, stS+arlen*ctS], 'stemWidth', wid);
end



% ah = annotation('arrow',...
%         'Color', [0 0 0],...
%         'headStyle','cback1','HeadLength',50,'HeadWidth',50);
% set(ah,'parent',gca);
% set(ah,'position',[X(ii) Y(ii) X(ii+1)-X(ii) Y(ii+1)-Y(ii)]);
% 
% % draw line along m
% quiver3(0, 0, 0, 0, -M-1, 0, 'k');
% quiver3(0, 0, 0, 0, M+1, 0, 'k');
% quiver3(0, 0, 0, 0, 0, 1.1*maxA, 'k');
% 
% thetaS = pi/2-pi/2/3;
% theta = linspace(0,thetaS,101);
% ys = zeros(size(theta));
% xs = 0.3*maxA*cos(theta);
% zs = 0.3*maxA*sin(theta);
% plot3(xs,ys,zs,'k');
% 
% arx = -sin(thetaS);
% arz = cos(thetaS);
% quiver3(0.3*maxA*cos(thetaS), 0, 0.3*maxA*sin(thetaS), arx, 0, arz, 'k')