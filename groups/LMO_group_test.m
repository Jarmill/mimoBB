rng(60)

w.weights = [2 1 1];
w.groups  = {[1 2], [3 4], 5};

G = randn([5, 1]);
norm_type = 2;

a = LMO_1d(G, norm_type, w);