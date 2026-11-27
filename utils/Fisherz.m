function x = Fisherz(x)
% Original finite Fisher transform, including finite values for +/-1.
assert(isreal(x), 'Fisherz expects real input.');
x(x > 1) = 1;
x(x < -1) = -1;
x = atanh(x);
x(isinf(x) & x > 0) = atanh(1 - eps);
x(isinf(x) & x < 0) = atanh(-1 + eps);
end
