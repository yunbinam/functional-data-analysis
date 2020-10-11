clear
addpath('..')
R1 = mmread('R1_642.mtx');

[X, beta, eigV] = mysamples(R1, 100, 50, 9, .01);
y = X*beta + normrnd(0, 1, [size(X,1),1]);
A = [y X];

save('eigV.mat', 'eigV');
save('beta_642.mat', 'beta');
save('X_samples_642.mat', 'X');
save('A_samples_642.mat', 'A');
