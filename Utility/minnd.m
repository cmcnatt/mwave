function [val, inds] = minnd(mat)

[val, inds] = maxnd(-mat);
val = -val;
