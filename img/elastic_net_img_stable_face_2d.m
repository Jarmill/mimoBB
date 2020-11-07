%shows stable-optimal, stable-suboptimal, and unstable faces in 2d

%matrix
%A = [1, 1; 0, 1];
%b = [1.5; 1.2];
%rng(9);
%rng(10);
%rng(23)
%rng(12)
%rng(35)

%rng(25);
%rng(51);
rng(93);


%A = 5*(rand(2, 2)-0.5);
%b = 7*(rand(2, 1)-0.5);

A = randn(2, 2);
b = 2.2*randn(2, 1);

x_ols = A \ b;

%regularization
delta = 1e-3;
tau = 1;
%coord_lim = [-1.25*tau 1.25*tau];
coord_lim = [-1.75*tau 1.75*tau];
%coord_lim = [-2*tau 2*tau];

%contour matrix
K = A'*A;
[~, S, V] = svd(K);
cont_matrix = V * sqrt(pinv(S));
t = linspace(0, 2*pi, 201);
circ = [cos(t); sin(t)];
contour_s = cont_matrix  * circ;

%ball
e1 = [1; 0];
e2 = [0; 1];
corners = [e1 e2 -e1 -e2]*tau;
ball = [corners tau*e1];

%extended face
fpp = [coord_lim; (tau-coord_lim)];
fnp = [coord_lim; (tau+coord_lim)];
fpn = [coord_lim; (-tau+coord_lim)];
fnn = [coord_lim; (-tau-coord_lim)];

%affine points
xpp = affine_inner([ 1;  1], tau, A, b, 0, 1, delta, [1; 1]);
xpn = affine_inner([ 1; -1], tau, A, b, 0, 1, delta, [1; 1]);
xnp = affine_inner([-1;  1], tau, A, b, 0, 1, delta, [1; 1]);
xnn = affine_inner([-1; -1], tau, A, b, 0, 1, delta, [1; 1]);

%collate
x_affine = [corners xpp xpn xnp xnn];
fval = en_abs(A, x_affine, b, delta);
fols = en_abs(A, x_ols, b, delta);
radius = sqrt(2*(fval - fols));

stable_face = sum(abs(x_affine), 1) <= tau;

x_stable = x_affine(:, stable_face);
x_unstable = x_affine(:, ~stable_face);


[~, imin] = min(fval(stable_face));
x_en = x_stable(:, imin);


%Gradients
grad = A'*(A*x_affine - b) + delta*x_affine;
grad_norm = grad ./ sqrt(sum(grad.^2, 1));

grad_norm_stable = grad_norm(:, stable_face);
grad_norm_unstable = grad_norm(:, ~stable_face);

%plots and output
figure(89)
clf
hold on
%start plotting

%ball and faces
plot(ball(1, :), ball(2, :),  'Color', [0.5, 0.5, 0.5], 'LineWidth', 2)
plot(fpp(1, :), fpp(2, :), ':k')
plot(fnp(1, :), fnp(2, :), ':k')
plot(fpn(1, :), fpn(2, :), ':k')
plot(fnn(1, :), fnn(2, :), ':k')

%error contours
Ncolor = 100;
par = parula(Ncolor+1);
fmin = min(fval);
fmax = max(fval);
y = floor(((fval - fmin)/(fmax - fmin))*Ncolor) + 1;

%remove saturation on contours of unstable faces
%may be a little too fancy
colors_raw = par(y, :);
colors_hsv = rgb2hsv(colors_raw);
colors_hsv(~stable_face, 2) = tanh(colors_hsv(~stable_face, 2)/3);
colors_hsv(~stable_face, 3) = tanh(colors_hsv(~stable_face, 3)*5);
colors_back = hsv2rgb(colors_hsv);

for i = 1:length(radius)
    curr = radius(i)*contour_s + x_ols;
    %plot(curr(1, :), curr(2, :), 'color', par(y(i), :))
    plot(curr(1, :), curr(2, :), 'color', colors_back(i, :))
end

%ols point
scatter(x_ols(1), x_ols(2), 200, 'xr')

grad_scale = 0.25;
%grad_scale = 0.4;
quiver(x_stable(1, :), x_stable(2, :), ...
    -grad_scale*grad_norm_stable(1, :), -grad_scale*grad_norm_stable(2, :), ...
    'filled', 'AutoScale', 'Off', 'color', [0.1, 0.6, 0])

quiver(x_unstable(1, :), x_unstable(2, :), ...
    -grad_scale*grad_norm_unstable(1, :), -grad_scale*grad_norm_unstable(2, :), ...
    'Filled', 'AutoScale', 'Off', 'color', [0.1, 0.7, 0.9])

scatter(x_en(1), x_en(2), 200, 'ok')

xlim(coord_lim);
ylim(coord_lim);
axis square
%axis off
hold off


figure(90)
clf

select = [1; 5; 6; 8]; %vertex, stable-opt, unstable, stable-subopt
titles = {'Vertex', 'Stable-Optimal', 'Unstable', 'Stable-Suboptimal'};
extended_face = {[0 0; 0 0], fpp, fpn, fnn};
for ti = 1:length(select)
    subplot(1, 4, ti)
    hold on
    %ball
    plot(ball(1, :), ball(2, :),  'Color', 0.5*[1 1 1], 'LineWidth', 2)
    scatter(x_ols(1), x_ols(2), 200, 'xr')
    scatter(x_en(1), x_en(2), 200, '*k')

    i = select(ti);
    
    %face
    face_curr = extended_face{ti};
    plot(face_curr(1, :), face_curr(2, :), '--k')
    
    %gradient
    x_curr = x_affine(:, i);
    grad_curr = 0.4*grad_norm(:, i);
%     arrow3D([x_curr; 0], -[grad_curr; 0], 'k', 0.6, 0.05, 3);

    %contour
    curr = radius(i)*contour_s + x_ols;
    plot(curr(1, :), curr(2, :), 'color', colors_raw(i, :), 'LineWidth', 2)

    hold off
    xlim(coord_lim);
    ylim(coord_lim);
    title(titles{ti})
    axis square
        
end