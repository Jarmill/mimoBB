function [A, w, poles] = pole_disk_grid(N_radius, N_horizon, group)
%POLE_DISK_GRID 
%   Sample a set of poles on the unit disk, and find their impulse
%   response for future fitting. Uniformly gridded sampling. A higher
%   N_points means a smaller grid size.

%radius = 40;
N_poles = 2*N_radius + 1; %number of poles on diameter of circle

rho = 1;

%Sample the upper half of the unit disk
[poles_xx, poles_yy] = meshgrid(linspace(-rho, rho, N_poles), linspace(0, rho, N_radius+1));
poles = poles_xx + 1.0j*poles_yy;
poles_circ = poles(abs(poles) <= 1);
%poles_circ = reshape(poles_circ, [length(poles_circ), 1]);


[ A, w, poles] = pole_matrix_upper( poles_circ, N_horizon, group);

end

