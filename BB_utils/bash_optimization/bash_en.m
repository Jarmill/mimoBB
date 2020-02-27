function [S_bash, c_bash, N_bag, error_bash, AS_bash, K_int_bash, s] = bash_en(A, b, tau, delta, CS, c0, S0, S_bag, AS0, K_int0, sparse_vec)
%BASH of bag and bash
%finds an optimal loading of atoms in S_bag to minimize LSQ + L2 problem
%(EN objective). Performs affine bashing with interior and/or exterior
%steps to solve for c_bash. Output is a stable-optimal face of the simplex
%or L1-ball in weight space.
%
%Input:
%   A, b, tau, delta:   parameters of objective EN problem (LSQ + L2)
%   CS:         If the total atomic set is centrally symmetric (0: simplex,
%               1: L1-ball). CS=1 allows for negative values of coefficients.
%   c0:         Initial set of coefficients found in the previous
%               iteration. Already optimal over its active set.
%   S0:         Set of atoms in previous active set, with optimal loadings c0.
%   S_bag:      Newly added atoms to system (representation). Added during
%               the BAG phase (cluster and/or orthogonal). If lockout 
%               occurs and not all of these can be added, half of the atoms
%               will be added. Atoms in S_bag are ranked in priority order.
%   AS0:        Matrix A times set of atoms in active set, computed in
%               previous iteration
%   K_int0:     Interior kernel matrix that was the output of a previous
%               bash step
%   s:          Signs of c_bash_out before flips and truncation. Optional.
%   sparse_vec: Basis is sparse vectors. Simplifies things
%
%Output:
%   S_bash:     atoms of active set with nonzero coefficients
%   c_bash:     nonzero coefficients of active set (zero values have been
%               bashed away)
%   N_bag:      the number of atoms from the bag that were successfully added
%   error_bash: Function value for f(Sc) at the end of bashing
%   AS_bash:    A*S_bash, quantity for next iteration
%   K_int_bash: K_int indexed through nonzero coefficients, precomputed for
%               next iteration

%Let's get this party started.

%general setup
N_prev = size(S0, 2);
N_bag_0 = size(S_bag, 2);
N_bag = N_bag_0;
%ind_bag = N_prev + (1:N_bag);
ind_bag = 1:N_bag;
ind_all = 1:(N_prev + N_bag_0);
done = 0; %done with BASH step

USE_CHOL = 0;

% Some atomic sets may have a redundant dictionary, and may add a new
% component which can be expressed as a linear combination of items in the
% bag. I think that's ok and to be expected. The kernel matrix K_int will
% then be PSD rather than PD. Add a small regularization bump to
% compensate.

%parameters for the EN-type problem
AS_bag = A*S_bag;
if isempty(S0)
    Kx = [];        %cross term on kernel
    rhs_int0 = [];  %right hand side of interior expression
else
    Kx =  AS0'*AS_bag + delta*S0'*S_bag;
    rhs_int0 = AS0'*b;
end

%setup kernel parameters
Kbag = (AS_bag'*AS_bag) + delta*(S_bag'*S_bag);

while ~done
    %atoms that are currently in the bag, and the set of atoms under
    %consideration with their initial weights. Atoms newly added from the
    %bag have a weight of 0.
    c = [c0; zeros(length(ind_bag), 1)];
    
     %the interior matrices can be used as a reference for exterior steps
     %build them up through progress from previous iterations
     ind_total = [1:N_prev (N_prev + ind_bag)];

    %interior kernels and solutions
    if isempty(S0)
        Kx_curr = [];
    else
        Kx_curr = Kx(:, ind_bag);
    end
    rhs_int_curr = AS_bag(:, ind_bag)'*b;
    rhs_int = [rhs_int0; rhs_int_curr];
    
    
    Kbag_curr = Kbag(ind_bag, ind_bag);
    
    K_int = full([K_int0 Kx_curr; Kx_curr' Kbag_curr]);
    %K_int   = (AS'*AS) + delta*(S'*S); %is r*r, where card(r) = r


    %K_int is singular when the atoms are linearly dependent. Figure out a
    %way to regularize against this, or another scheme to take this into
    %account. Entries of K_int are usually large, so this may be used to
    %adjust epsilon?
    %epsilon = 1e-6;    
    if sparse_vec
        %built in orthogonal set of atoms, makes life easy.
        epsilon = 0;
    else
        epsilon = 1e-10*max(diag(K_int)); %adaptive regularization?
    end
        
    %with cholesky
    if USE_CHOL
        %with cholesky
        c_bash_out = bash_en_inner_chol(c, tau,  K_int + eye(size(K_int))*epsilon, rhs_int, CS);
        done = 1;
    else
        %without cholesky    
        [c_bash_out, lockout, lockout_dims, on_boundary] = bash_en_inner(c, tau, K_int + eye(size(K_int))*epsilon, rhs_int, CS);    
        %end of processing
        if lockout
            %standard AIMD lockout prevention
            %N_bag = ceil(N_bag/2);
            %ind_bag = ind_bag(1:N_bag);

            %New lockout prevention, kill the locked dimensions
            ind_bag(lockout_dims - N_prev) = [];
            N_bag = length(ind_bag);
            if (N_bag_0 == 1) && lockout
                break
            end
        else
            done = 1;
        end
    end

end

%drop atoms with zero weight from active set
iactive = find(c_bash_out);
iactive_atom = ind_total(iactive);

%signs of everything
%make sure to flip the signs of S and AS if c is negative
s = zeros(size(ind_all));
s_active = sign(c_bash_out(iactive));
%full set of signs sent to atom_svd_bash to revise and knock out dimensions
s(iactive_atom) = s_active;

ind_old_in = iactive_atom(iactive_atom <= N_prev);
ind_new_in = iactive_atom(iactive_atom > N_prev)-N_prev;

%figure out nonzero atoms corresponding to entries in S and AS
S_old_in = S0(:, ind_old_in);
AS_old_in = AS0(:, ind_old_in);

S_new_in = S_bag(:, ind_new_in);
AS_new_in = AS_bag(:, ind_new_in);

N_bag = length(ind_new_in);

%collect all nonzero atoms for next iteration
S_bash = [S_old_in S_new_in];
AS_bash = [AS_old_in AS_new_in];
K_int_bash = K_int(iactive, iactive); %this cheese-clothing may be inefficient
c_bash = c_bash_out(iactive);
%flip the signs of atoms and their influence AS. Saves bagging on the
%boundary, and will be used for the cholesky black magic in BASH
if any(s_active == -1)
    %
    %bsxfun turns things into sparse vectors. Do I want that?
    %
    %for flipped values, where D = diag(s_active)
    %Sf = S D
    %ASf = AS D
    S_bash  = full(bsxfun(@times,S_bash,s_active'));
    AS_bash = full(bsxfun(@times,AS_bash,s_active'));
    %Kf = D K D
    %may be a more efficient way to handle this sign flipping
    K_int_bash = full(bsxfun(@times, s_active, bsxfun(@times, K_int_bash, s_active')));
    %all coefficients become positive now. Thank you central symmetry.
    c_bash = abs(c_bash);
end

error_bash = 0.5*(norm(AS_bash*c_bash - b)^2 + delta*norm(S_bash*c_bash)^2);
%error_bash = norm(AS*c_bash - b)^2 + delta*norm(S*c_bash)^2;