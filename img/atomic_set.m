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
scatter(target_x, target_y, 200,  'xk', 'linewidth', 3)
axis off
hold off

