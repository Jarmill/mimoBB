load('ring_sys_plt.mat')

RandRounds = opt.RandomRounds;

cost = out.cost_list;
card = out.card_list;
order = out.order_list;

figure(20)
clf
hold on
yyaxis left
%plot(card)
plot(order)
ylabel('System Order')
ylim([0, 24]);

FS = 18;
text(10, 22, 'Randomize', 'FontSize', FS)
text(22, 22, 'Reweight', 'FontSize', FS)


yyaxis right
plot(cost)
plot((RandRounds)*[1,1], ylim, 'k:')
ylabel('Fitting Error')
hold off

xlabel('Active-Set Iteration')

title('MIMO Identification and Sparsification')

figure(21)
clf
iopzplot(out.sys_out);