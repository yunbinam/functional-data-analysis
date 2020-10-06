clear
addpath('brain/regression test file')
R1 = mmread('R1.mtx');

[X, beta] = mysamples(R1, 100, 50, 9, .01);
y = X*beta + normrnd(0, 1, [size(X,1),1]);
A = [y X];

save('beta.mat', 'beta');
save('X_samples.mat', 'X');
save('A_samples.mat', 'A');
