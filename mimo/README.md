This directory contains functions and helpers that implement MIMO system identification. The main function to implement system identification is `atomic_LTI_iteration`. The options data structure is located in `mimoAtomOptions`. 

MIMO system identification is performed in two phases: Randomize and Reweight. 

The Randomize phase adds and selects new poles to the system, where the set of poles are sampled at random from the unit disk (`createAtoms` which calls `uniform_over_ring_sector`). The randomize phase uses an L1/Linf penalty where the group weights are scaled by the pole locations (`get_Scales`).

The Reweight phase sparsifies the selection of poles by performing a Reweighted Heuristic and scaling the atomic ball. No new atoms are added in the Reweight phase.

Each Randomize or Reweight iteration in `atomic_LTI_iteration` involves a call to `atomic_LTI_BB`. The function `atomic_LTI_BB` creates data operators for the Fully Corrective Frank-Wolfe (FCFW) code in `BB_operator'. `atomic_LTI_BB` then formats the solved coefficients/problem from the FCFW routine and recovers the solved system.

