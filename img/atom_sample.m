opt = mimoAtomOptions;
opt.NumAtoms = 10;
Ns = 50;

[ha,p, scales, groups, f, L] = createAtoms(Ns,opt);

figure(1)
clf
plot(ha)