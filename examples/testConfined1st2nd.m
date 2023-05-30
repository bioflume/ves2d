clear; clc;

Ns = [16; 32; 64; 128];
ave_iters = zeros(size(Ns));
max_iters = zeros(size(Ns));
for in = 1 : numel(Ns)
  [ave_iters(in), max_iters(in)] = confined_2ndKind(Ns(in),0);
%   [ave_iters(in), max_iters(in)] = confined_1stKind(Ns(in));
end