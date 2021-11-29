A modification of Marina Vinyes' fcgan code for fully corrective frank wolfe

For the atom system.

Hopefully this code avoids the regeneration issue in group sparsity that has been plaguing the BB implementation.




fcgan always remains on the face of a simplex, and solves an active set qp
routine in order to find optimal weightings within the current set of atoms.

Interior points are represented in an atomic set as (c+) a + (c-) (-a)
for positive and negative weights c+ and c-.

The hessian matrix in the QP solver is therefore rank deficient.

There is an issue where fcgan will select the same atom multiple times. 
This occurs when there are multiple atoms that return the same polar norm
value in the LMO. For the L1 case, this will happen when multiple coordinates
have the same absolute value of the gradient. The `max' call in the LMO 
(LMO_1d by my implementation) selects the atom with the lowest index, as is
standard in MATLAB. I don't know how to get rid of this, outside of performing
'unique' calls and selecting atoms that are outside the current set (which
is slow).