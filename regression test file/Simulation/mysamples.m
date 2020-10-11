function [X, beta, V] = mysamples(R1, n, m, v, noise) 
%
% Function: mysamples(R1,n,m)
% simulation sample generator
% outputs X: sample matrix, beta: smooth coefficients
% 
% Required arguments
%%% R1: stiff matrix
%%% n: the number of eigen functions
%%% m: the number of samples
%%% v: variance of random coefficients
%
% Optional arguments
% noise: random noise on each vertice
%

if (nargin == 4)
    noise = 0;
end

% generate smooth functions on the manifold
[V, D] = eigs(R1, n, 'smallestabs');

beta = V(:,round(n/2));
X = zeros(1, length(R1));

for i = 1:m
    r = randi([1 n],[1, 10]);
    e = V(:,r);
    e = e';
    c = normrnd(0,v,[1, 10]);
    x = c*e + normrnd(0,noise,[1,length(R1)]);
    
    X(i,:)=x;
end

end

