%demonstration of reweighting heuristic

%very important
%correct sign combo is [1,1]

% rng(1);
% rng(10);
% rng(23)
% rng(50);
% rng(60);
rng(48);
% 
% A = randn(2, 2);
% A(:, [1 2])=A(:, [2 1]);
% b = 2.2*randn(2, 1);
% b = 1.5*randn(2, 1);

% x_ols = A \ b;

% x_ols = [0.5; 1.4];
% x_ols = [0.5; 1.05];
% x_ols = [0.3; 1.4];
x_ols = [1.4; 0.3];
% x_ols = [0.4; 1.5];
A = randn(2, 2);
A(:, [1 2])=A(:, [2 1]);
b = A*x_ols;

%regularization
delta = 1e-3;
tau = 1;
%coord_lim = [-1.25*tau 1.25*tau];
coord_lim = [0 1.75*tau];
%coord_lim = [-2*tau 2*tau];

%contour matrix
K = A'*A;
[~, S, V] = svd(K);
cont_matrix = V * sqrt(pinv(S));
t = linspace(0, 2*pi, 201);
circ = [cos(t); sin(t)];
contour_s = cont_matrix  * circ;


%extended face
[x_en, lam_en] = aff_pos2(A, b, delta, [1; 1], tau);
w = 1./(abs(x_en) + 1e-4);
w = tau * w / (w' * x_en);
% w = tau./(2 * x_en);
[x_r, lam_r] =  aff_pos2(A, b, delta, w, tau);


w2 = 1./(abs(x_r) + 1e-4);
w2 = tau * w2 / (w2' * x_r);
% w = tau./(2 * x_en);
[x_r2, lam_r2] =  aff_pos2(A, b, delta, w2, tau);

grad_en = A'*(A*x_en - b) + delta*x_en;
grad_r = A'*(A*x_r - b) + delta*x_r;
grad_r2 = A'*(A*x_r2 - b) + delta*x_r2;



fols = en_abs(A, x_ols, b, delta);
fen = en_abs(A, x_en, b, delta);
fr   = en_abs(A, x_r, b, delta);
fr2   = en_abs(A, x_r2, b, delta);


rad_en = sqrt(2*(fen- fols));
rad_r = sqrt(2*(fr- fols));
rad_r2 = sqrt(2*(fr2- fols));
% radius = [rad_en; rad_r];
radius = [rad_en; rad_r; rad_r2];

figure(91)
clf
hold on

%ball
e1 = [1; 0];
e2 = [0; 1];
corners = [e1 e2 -e1 -e2]*tau;
plot([tau; 0], [0; tau], 'k', 'LineWidth', 2)
plot(1/w(1)*[tau; 0], 1/w(2)*[0; tau], '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2)
plot(1/w2(1)*[tau; 0], 1/w2(2)*[0; tau], ':', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 2)
% ball = [corners tau*e1];
% plot(ball(1, :), ball(2, :),  'Color', [0.5, 0.5, 0.5], 'LineWidth', 2)
% 
% %faces
% fpp = [coord_lim; (tau-coord_lim)];
% fnp = [coord_lim; (tau+coord_lim)];
% fpn = [coord_lim; (-tau+coord_lim)];
% fnn = [coord_lim; (-tau-coord_lim)];
% % 
% plot(fpp(1, :), fpp(2, :), ':k')
% plot(fnp(1, :), fnp(2, :), ':k')
% plot(fpn(1, :), fpn(2, :), ':k')
% plot(fnn(1, :), fnn(2, :), ':k')


% plot(1/w(1)*fpp(1, :), 1/w(2)*fpp(2, :), ':k')
% plot(1/w(1)*fnp(1, :), 1/w(2)*fnp(2, :), ':k')
% plot(1/w(1)*fpn(1, :), 1/w(2)*fpn(2, :), ':k')
% plot(1/w(1)*fnn(1, :), 1/w(2)*fnn(2, :), ':k')

%optimal points
scatter(x_ols(1), x_ols(2), 200, 'xr')
scatter(x_en(1), x_en(2), 200, 'ok')
scatter(x_r(1), x_r(2), 200, 'sk')

scatter(x_r2(1), x_r2(2), 200, '^k')

for i = 1:length(radius)
    curr = radius(i)*contour_s + x_ols;
    %plot(curr(1, :), curr(2, :), 'color', par(y(i), :))
    plot(curr(1, :), curr(2, :))
end
xlim([0; 1.75]);
ylim([0; 1]);

text(x_en(1)-0.05, x_en(2)-0.05, 'x0', 'FontSize', 14, 'Interpreter', 'Latex')
text(x_r(1)+0.02, x_r(2)+0.05, 'x1', 'FontSize', 14, 'Interpreter', 'Latex')
text(x_r2(1)+0.02, x_r2(2)+0.05, 'x2', 'FontSize', 14, 'Interpreter', 'Latex')
% axis square
hold off

axis off

% [x_en x_r x_r2]

function [x, lam] =  aff_pos2(A, b, delta, w, tau)
    %solve elastic net problem on 2d ball [1, 1]
    K = A'*A + delta*eye(length(w));
    ref = A'*b;
    
    KKT = [K w; w' 0];
    KKTref = [ref; tau];
    
    KKTout = KKT \ KKTref;
    
    x = KKTout(1:2);
    lam = KKTout(3);
end
