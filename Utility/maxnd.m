function [val, inds] = maxnd(mat)

N = ndims(mat);
Inds = cell(N, 1);

[valsn, indsn] = max(mat);
Inds{1} = indsn;

for n = 2:N
    [valsn, indsn] = max(valsn);
    Inds{n} = indsn;
end

val = valsn;
inds = zeros(1, N);

inds(N) = Inds{N};

for n = (N-1):-1:1
    matn = squeeze(Inds{n});
    cinds = num2cell(inds((n+1):N));
    inds(n) = matn(cinds{:});
end