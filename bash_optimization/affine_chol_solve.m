function x_out = affine_chol_solve(B, p, K_full, L_full, rhs_full)
%Cholesky solver for use in affine steps. Solves for x_full and x_aff in a
%very efficient manner. When K = L'L, solves KB x = r for x.
%
%Input:
%   B:      Affine subspace directions, so c = By+p
%   p:      Affine subspace base point
%   K_full:  Interior kernel matrix among selected atoms
%   L_full:  Cholesky factorization of K_full (interior kernel matrix)
%   rhs_aff_partial:  answer vector of exterior = (rhs_full - Kp)
%
%not exactly cholesky, but close enough.
%who knows, maybe it is fast.


rhs_aff = B'*(rhs_full - K_full*p);

%raw, used for reference
%K_aff = B'*K_full*B;
%y_aff_raw = K_aff \ rhs_aff;

%actual processign
%[Q_aff, L_aff] = qr(L_full*B, 0);
[~, L_aff] = qr(L_full*B, 0);

%L_aff = qr(L_full*Be, 0);
%norm(Q_aff*L_aff - L_full*Be)


y_aff_partial = L_aff' \ rhs_aff;
y_aff = L_aff \ y_aff_partial;

%norm_diff = norm(y_aff_raw - y_aff);

%return output
%may be sensitive to content in left nullspace
x_out = B*y_aff + p;
end