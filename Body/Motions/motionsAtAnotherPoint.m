function [motOut] = motionsAtAnotherPoint(motIn, cg, point)

[Nt, Nb, dof] = size(motIn);

if dof < 6
    error('motionsAtAnotherPoint requires at least 6 DOF. 1,2,3 are surge, sway, heave; 4,5,6 are roll, pitch, yaw');
elseif dof > 6
    warning('motionsAtAnotherPoint assumes that DOF 1,2,3 are translation modes and DOF 4,5,6 are rotational modes');
end

dx = point(1) - cg(1);
dy = point(2) - cg(2); 
dz = point(3) - cg(3);

R = [0      dz     -dy;
    -dz     0       dx;
     dy    -dx      0];

motOut = complex(zeros(size(motIn)));

for m = 1:Nt
    for n = 1:Nb
        motOut(m, n, 1:3) = squeeze(motIn(m, n, 1:3)) + (R*squeeze(motIn(m, n, 4:6)));
        motOut(m, n, 4:end) = motIn(m, n, 4:end);
    end
end

