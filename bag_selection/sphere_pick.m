function [x, y, z] = sphere_pick(N)
%SPHERE_PICK picks N points on surface of sphere with radius 1
%outputs: x, y, z coordinates of sphere

%uniform random variables
z = 2*rand(1, N)-1;
v = rand(1, N);

theta = 2*pi*v;

x = sqrt(1 - z.^2) .* cos(theta);
y = sqrt(1 - z.^2) .* sin(theta);
