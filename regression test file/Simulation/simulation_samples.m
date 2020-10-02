clear
addpath('brain/regression test file/Simulation')
R1 = mmread('R1.mtx');

[X, beta] = mysamples(R1, 100, 50, 9, .01);

save('beta.mat', 'beta');
save('X_samples.mat', 'X');
