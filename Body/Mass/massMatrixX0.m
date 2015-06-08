function [M0] = massMatrixX0(M, x0)
% computes the 6 dof mass matrix, M0, about some point x0, where M is the 
% orginal 6 dof mass matrix computed about its cg, and x0 is a 3D vector 
% given in the coordinates of M

M0 = M;

mass = M(1,1);

dx = x0(1);
dy = x0(2);
dz = x0(3);

dmat = [0 -dz dy; dz 0 -dx; -dy dx 0];

% parallel axis theorem for modifying moment of interia
dmat2 = dmat*dmat;

M0(4:6, 4:6) = M0(4:6, 4:6) - mass*dmat2;

% moments due to cg torque
M0(1:3, 4:6) = mass*dmat;
M0(4:6, 1:3) = -mass*dmat;



