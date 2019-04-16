% Size of the random variable
M = 10;
% Number of time samples
N = 1000;
% Number of monte carlo trials
T = 10000;

% sample random vector to make matrix
z = randn(M,1);
A = toeplitz(z);

% Calculate eigenvalues
[V,D] = eig(A);
val = diag(D);

% Shift matrix if negative eigenvalues
lamb_min = min(val);
th = z;
if lamb_min < 0
  % add some uniform random noise to make it strictly positive definite
  th(1) = th(1) + abs(lamb_min) + rand();
end

% Construct correlation matrix
R = toeplitz(th);

e = zeros(T,M, M);
for t = 1:T
  % Generate random samples
  sqrtR = sqrtm(R);
  Y = rand(M,N);
  X = sqrtR*Y;

  Rhat = X*X'/N;
  thhat = Rhat(:,1);
  e(t, :, :) = (th - thhat)*(th - thhat)';
end

sigma2e = sum(e)/T;
sigma2e = reshape(sigma2e, M, M);

% compute CRB

KR = kron(R,R);
jac = [];
I10 = eye(10);

for i = 1:M
  e_i = I10(:,i);
  A_i = toeplitz(e_i);
  jac = [jac A_i(:)];
end

J = (N/2)*jac'*KR*jac;
C = inv(J);
