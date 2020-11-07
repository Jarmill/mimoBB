N = 8;

rng(30, 'twister')

th = (1:N)*(2*pi)/N;

b = 0.3 * (2*rand(1, N)-1);
th_off = th + b;

th_off2 = [th_off th_off(1)];


r = 1.5;
alpha = 0.7;
ind = 2;
target_x = r*(alpha*cos(th_off(ind)) + (1-alpha)*cos(th_off(ind+1)));
target_y = r*(alpha*sin(th_off(ind)) + (1-alpha)*sin(th_off(ind+1)));



figure(1)
clf
hold on
plot(r*cos(th_off2), r*sin(th_off2))
plot(cos(th_off2), sin(th_off2), '--k')
scatter(cos(th_off2), sin(th_off2), 300, '.k')
scatter(target_x, target_y, 400,  'xk', 'linewidth', 3)
axis off
hold off

figure(2)
clf
grad = 0.9*randn(2, 1);
grad = grad/norm(grad);
lmo = grad'*[cos(th_off); sin(th_off)];
[m, i] = max(lmo);
hold on
plot(cos(th_off2), sin(th_off2), '--k')
quiver(0,0,grad(1), grad(2), 'LineWidth', 5, 'maxheadsize', 0.6)
scatter(cos(th_off2), sin(th_off2), 300, '.k')
scatter(cos(th_off(i)), sin(th_off(i)), 300, 'ok')
text(-0.45, -0.4, '$\nabla f(x)$', 'interpreter', 'latex', 'fontsize', 28)
axis off
hold off