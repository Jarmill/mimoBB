function N_next = bag_size_update(N_survived, N_dropped, int_increase, ext_increase, on_boundary, on_boundary_new, N_max)
%Feedback control of the bag size, number of atoms to be added in BAG
%phase. Inspired by AIMD (additive increase/multiplicative decrease) of
%TCP.
%
%Input:
%   N_survived:     Number of atoms from bag that survived BASH
%   N_dropped:      Number of atoms from previous set that survived BASH
%   int_increase:   Bag size increase if c on interior
%   ext_increase:   Bag size increase if c on exterior
%   on_boundary:    Whether old set of coefficients on exterior
%   on_boundary_new:Whether new set of coefficients on exterior
%   N_max:          Maximum allowable bag size

if nargin < 6
    N_max = Inf;
end

if on_boundary ~= on_boundary_new
    %switch from interior to exterior
    if on_boundary_new
        N_next = max(ext_increase, 1);
    else
        N_next = int_increase;
    end
else    
    %maintains partition (stays in interior or on exterior)
    if on_boundary_new
        N_next = N_survived + ext_increase;
    else
        N_next = N_survived + int_increase;
    end
    
    %wind down bagging if too many atoms are being dropped
    %need a better strategy
    if N_dropped > N_survived
        N_next = ceil(N_next/2);
    end
    
    if N_next > N_max
        N_next = N_max;
    end
end

end