rng(89)
%nInstances = 500;
%nVars = 50;
nInstances = 1000;
nVars = 100;

%tau = 2;
% nInstances = 4000;
% nVars = 500;
X = randn(nInstances,nVars);
z = complex(randn(nVars,1),randn(nVars,1)).*(rand(nVars,1) > .5);
y = X*z;

%tau = 100;
tau = 2;
delta = 0;

[x_bb, S_bb, c_bb] = BB_1d(X, y, tau, delta, 1);
%[x_bb, S_bb, c_bb] = BB_1d(X, y, tau, delta, Inf);
% 
% 
% 
% grad_bb = X'*(X*x_bb - y) + delta*x_bb;
% 
% figure(4)
% clf
% 
% N = size(x_bb, 1);
% ind_range = 1:N;
% 
% %plot3([1, N], [0, 0], [0, 0], 'k', 'LineWidth', 1)
% th = linspace(0, 2*pi, 51);
% circ = [cos(th); sin(th)];
% %scatter3(, real(grad_bb), imag(grad_bb))
% subplot(2, 1, 1)
% hold on
% scatter3(ind_range(x_bb ~= 0), real(grad_bb(x_bb ~= 0)), imag(grad_bb(x_bb ~= 0)), 100, '.k')
% scatter3(ind_range(x_bb == 0), real(grad_bb(x_bb == 0)), imag(grad_bb(x_bb == 0)), 100, [21 200 225]/255.0, '.')
% 
% for i = ind_range    
%     grad_curr = grad_bb(i);
%     r = abs(grad_curr);
%     if x_bb(i) ~= 0
%         plot3([i, i], [0, real(grad_curr)], [0, imag(grad_curr)], 'k')
%         plot3(ones(size(th))*i, r*circ(1, :), r*circ(2, :), ':k')
%     else
%         plot3([i, i], [0, real(grad_curr)], [0, imag(grad_curr)], 'Color', [21 200 225]/255.0)
%     end
% end
% title('Complex Gradient at convergence')
% xlabel('$j$', 'Interpreter', 'latex')
% ylabel('Re $\nabla f(x)$', 'Interpreter', 'latex')
% zlabel('Im $\nabla f(x)$', 'Interpreter', 'latex')
% hold off
% view(3)
% %axes.SortMethod='ChildOrder';
% 
% subplot(2, 1, 2)
% x_prev = x_bb * exp(1j*pi/6);
% hold on
% x_prev_active = find(x_prev);
% for i = ind_range
%     x_curr = x_prev(i);
%     plot3([i, i], [0, real(x_curr)], [0, imag(x_curr)], 'color', [0    0.4470    0.7410])    
% end
% scatter3(x_prev_active, real(x_prev(x_prev_active)), imag(x_prev(x_prev_active)), 100, [0    0.4470    0.7410],  '.');
% 
% 
% 
% 
% x_active = find(x_bb);
% for i = ind_range
%     x_curr = x_bb(i);
%     plot3([i, i], [0, real(x_curr)], [0, imag(x_curr)], 'color', [0.8500    0.3250    0.0980])    
% end
% scatter3(x_active, real(x_bb(x_active)), imag(x_bb(x_active)), 100, [0.8500    0.3250    0.0980],  '.');
% 
% 
% hold off
% 
% hold on
% N_atoms = size(S_bb, 2);
% for k = 1:N_atoms
%     c_curr = c_bb(k);
%     i_curr = find(S_bb(:, k));
%     for i = i_curr
%         S_curr = S_bb(i, k);
%         plot3([i, i], tau*c_curr*[0, real(S_curr)], tau*c_curr*[0, imag(S_curr)], 'k')
%         plot3([i, i], tau*[0, real(S_curr)], tau*[0, imag(S_curr)], ':k')
%         plot3(ones(size(th))*i, tau*circ(1, :), tau*circ(2, :), ':k')
%     end
% end
% hold off
    %plot(abs(grad_bb))

