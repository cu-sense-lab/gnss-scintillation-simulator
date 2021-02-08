function nfft = nicefftnum(n)
% Finds the smallest even number `nfft` that is greater than or equal to `n`
%  that is divisible by composite odd numbers
testvalues = [1, 3, 5, 9, 15, 21, 25, 27, 33, 35, 39];
nt = 2.^max(1, ceil(log2(n ./ testvalues))) .* testvalues;
nfft = min(nt);
return
