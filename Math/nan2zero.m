function [Mout] = nan2zero(Min)

Mout = Min;

inds = isnan(real(Min));
Mout(inds) = 0 + 1i*imag(Mout(inds));

inds = isnan(imag(Min));
Mout(inds) = real(Mout(inds)) + 1i*0;

end