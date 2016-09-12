function [yOut, f, S, fcut] = fftFilter(x, y, varargin)

L = length(y);
dx = x(2) - x(1);
N = floor(L/2);

% fourier transform
Y = fft(y);

% frequencies
f = (0:N)/dx/L;
% spectrum
S = 2/L*Y(1:N+1);

% sample frequency and cutoff frequency
fsamp = 1/dx;
if isempty(varargin)
    fcut = fsamp/8;
else
    fcut = varargin{1};
end

[~, indf] = min(abs(f - fcut));

Y2 = zeros(size(Y));
Y2(1:indf) = Y(1:indf);
% add the complex conjugate portions for the ifft
inds = length(Y2):-1:(length(Y2)-indf + 2);
Y2(inds) = Y(inds);

yOut = ifft(Y2);

