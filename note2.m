% -------- Formulating numerical solution of Orr-Sommerfeld Eqn ----------
% Solving the generalized EVP with Re_c = 5772.22 and alpha=1.02056 should 
% yield an eigenvalue which has c_i approx 0. This corresponds to the most 
% unstable mode, the onset of instability that is.
%%
c = [];
Re = 5772.22;
alpha = 1.02056;
for m=12:12
    N = 2^m;
    gridpts = 2^(m-1); % # grid points
    h = 2/(gridpts-1);  
    U = zeros([gridpts 1]);
    for i = 1:gridpts-1
        U(i)= 1-(-1 + h*(i-1))^2;
    end
    A = zeros(N);
    B = zeros(N);
    b = (-1)^0.5 / (alpha*Re);
    % phi' = 0 BC
    A(1,2) = -3;
    A(1,4) = 4;
    A(1,6) = -1;
    % phi' = 0 BC
    A(2,2) = 1;
    for j = 2:N/2 - 1
        k = 2*j-1;
        B(k,k) = 1;
        A(k,k-2) = b/h^2 ; %psi_j-1
        A(k,k) = U(j) - (2/h^2 + alpha^2)*b ; %psi_j
        A(k,k+1) = 2 ; %phi_j
        A(k,k+2) = b/h^2; %psi_j+1
        k = k+1;
        A(k,k-2) = -1/h^2;
        A(k,k-1) = 1;
        A(k,k) = 2/h^2 + alpha^2; 
        A(k,k+2) = -1/h^2;
    end
    % phi' = 0 BC
    A(N-1,N-4) = 1;
    A(N-1,N-2) = -4;
    A(N-1,N) = 3;
    % phi = 0 BC
    A(N,N) = 1;
    [eigvec,v] = eig(A,B); % v generalized eigenvalues as entries on a diagonal matrix 
    vsorted = sort(diag(v));
    vsorted(1) % Least stable mode
    gridpts; % # grid points
    c = [c;abs(imag(vsorted(1)))];
end
v = diag(v);

