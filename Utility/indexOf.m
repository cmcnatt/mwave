function [ind] = indexOf(vect, val)

[~, ind] = min(abs(vect - val));