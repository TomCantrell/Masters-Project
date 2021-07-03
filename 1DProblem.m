
% ===============================================================
% ------------    Linear 1D Eigenvalue Problem    ---------------
% ===============================================================

% Solving the linear one-dimensional eigenvalue problem 

k = 5; % First k Eigenvalues 
p = 11;
eig_exact = zeros([k 1]);
M = zeros([k p-2]);
H = [];
for i = 1:k
    eig_exact(i) = i^2 * pi^2; 
end

for m=3:10
  N = 2^m
  h = 1/(N-1); % N nodes => 1/N-1 step size for a domain of size 1
  H = [H;h];
  A = zeros(N);
  B = zeros(N);
  A(1,1) = -1.0/h^2;  % f_1=0
  for i = 2:N-1  % only one loop is needed
    B(i,i) = 1;
    A(i,i)= 2/h^2;
    A(i,i+1) = -1/h^2;
    A(i,i-1) = -1/h^2;
  end
  A(N,N) = -1.0/h^2;  % f_N=0
  [V,W] = eig(A,B); % Eigenvalues not in order
  v = sort(diag(W));
  % create vector of approx first k eigenvalues
  eig_approx = v(3:k+2)
  error = abs(eig_approx - eig_exact)
  M(:,m-2:m-2) = error;
end
for j=1:N
   if W(j,j)==eig_approx(3)
       eig_fn = V(1:N,j); 
   end
end
z = linspace(0,1,N);
plot(z,eig_fn,"color","k","LineWidth",1.0)
grid on 
xlabel("z")
ylabel("f")
xlim([-0.1 1.1])
ylim([-1.1 1.1])

%% create plot
gradient = zeros([k 1]);
for j = 1:k
    loge = log(M(j:j,:));
    logH = log(H);
    m = (loge(length(loge))-loge(1)) / (logH(length(logH))-logH(1))
    gradient(j) = m;
    p = plot(logH,loge);
    hold on
end
xlabel("log(h)")
ylabel("log(|\lambda^2 - i^2\pi^2|)")
%title("Log plot of Error of Eigenvalues vs step size h.")
legend("\lambda_1","\lambda_2","\lambda_3","\lambda_4","\lambda_5","Location","northwest") 
grid on 







