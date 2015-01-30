function [am] = cirWaveCoefsFromFunc(a, M)

N = length(a);

fori = 2*pi/N*fft(a);

am = ones(2*M+1,1);

am(M+1) = fori(1);
am(M+2:2*M+1) = fori(2:M+1);
am(M:-1:1) = fori(N:-1:N-M+1);

m = (-M:M).';
am = (-1i).^m.*am;

am = am.';
end