function omega = recover_via_exprit(X, Mh)

[U Lambda] = eig(X); 

N = size(X, 1);
% Estimate the number of harmonics
I = N:-1:N-Mh+1;


Y = U(:, I)*sqrt(Lambda(I, I));
GY = Y(1:end-1, :);
FY = Y(2:end, :);
LambdaS = GY\FY;
[Q LambdaSS] = eig(LambdaS);
omega = angle(diag(LambdaSS))/(2*pi);
omega = sort(mod(omega,1));

