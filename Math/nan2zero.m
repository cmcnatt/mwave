function [Mout] = nan2zero(Min)

Mout = Min;
inds = isnan(Min);
Mout(inds) = 0;

end