opt = mimoAtomOptions;
opt.NumAtoms = 130;
opt.IncludeConstant = 0;

rng(25)

[h,p_out,scales, groups, f, L] = createAtoms(60,opt);

opt2 = mimoAtomOptions;
opt2.NumAtoms = 50;
opt2.r1 = 0.7;
opt2.r2 = 0.9;
opt2.phi2 = pi/3;
opt2.IncludeConstant = 0;

[h2,p_out_2,scales2, groups2, f2, L2] = createAtoms(60,opt2);
a2 = angle(p_out_2);

p2 = p_out_2((a2 < pi/2) & (a2 > -pi/2));
% p2 = p_out_2((
th =  linspace(0, 2*pi, 200);


figure(1)
clf
subplot(1,2,1)
hold on
scatter(real(p_out), imag(p_out), 100, 'xk')
plot(cos(th), sin(th), 'k', 'LineWidth', 2)
xlim([-1,1])
ylim([-1,1])
axis square
axis off
title('Unit Disk','Fontsize', 16)

subplot(1,2,2)
hold on
scatter(real(p2), imag(p2), 100, 'xk')
plot(cos(th), sin(th), 'k', 'LineWidth', 2)
xlim([-1,1])
ylim([-1,1])
axis square
axis off
title('Sector Bound', 'Fontsize', 16)