% Size of the random variable
M = 10;
% Number of time samples
N = 80;
% Number of monte carlo trials
T = 1000;

% Make jacobian matrix
jac = [];
I10 = eye(10);

for i = 1:M
  e_i = I10(:,i);
  A_i = toeplitz(e_i);
  jac = [jac A_i(:)];
end

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
options = optimset('Display', 'off');

e = zeros(T,M, M);
eS = zeros(T,M, M);
for t = 1:T
  % Generate random samples
  sqrtR = sqrtm(R);
  Y = randn(M,N);
  X = sqrtR*Y;

  S = X*X'/N;
  f = @(theta) regression(theta,S,jac);
  th0 = zeros(M,1);
  thhat = fsolve(f, th, options);
  thhatS = S(:,1);
  e(t,:,:) = (th - thhat)*(th - thhat)';
  eS(t,:,:) = (th - thhatS)*(th - thhatS)';
end

C = sum(e)/T;
C = reshape(C, M, M);
CS = sum(eS)/T;
CS = reshape(CS, M, M);

% compute CRB

KR = kron(R,R);

J = (N/2)*jac'*KR*jac;
CRB = inv(J);
idx = 0:M-1;
plot(idx,diag(CRB),'g', idx, diag(C), 'b', idx, diag(CS), 'r');
xlabel('Element Number');
ylabel('Error Variance');
legend('CRB', 'C', 'C_S');

function F = regression(theta,S,jac)
  R = toeplitz(theta);
  Rinv = inv(R);
  A = Rinv*(R-S)*Rinv;
  F = jac'*A(:);
end
