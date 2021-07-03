
% Implementation of plane Couette-Poiseuille flow

% Using spectral methods, only changes are to the base profile

%% Formulate the generalized eigenvalue problem 
Re = 6000; % Reynolds number
alpha = 1.097; % Wavenumber
u_w=0.0; % wall speed
N=54; % N+1 Chebyshev polynomials
% Working out E matrix 
A = zeros(N);
B = zeros(N);
C = zeros(N);
for i=1:N
    if i==1
        C(i,i)=2;
    else 
        C(i,i)=1;
    end
    if i<=N-2
        C(i,i+2)=-1;  
    end
    A(i,i+1)=2*i;
end
A; 
C; 
E = C\A;
E(N+1,N)=0; % Add another row to make the matrix square

% Constructing the matrix which is used to determine the matrix formulation
% of the cross terms in the Orr-Sommerfeld equation U*phi, U*phi"
Q = zeros(N+1); % has same dimensions as E
for j=1:N-3
    if j<N
        Q(j,j+2) = -0.375;
    end
    % Diagonal entries
    if j-1 == 1
        Q(j,j)=0.375;
    else
        Q(j,j)=0.75;
    end        
    % Other off-diagonal entries
    if j>2
        Q(j,j+1)=0.5*u_w;
        if j==3
            Q(j,j-2) = -0.75;
        else
            Q(j,j-2)=-0.375;
        end
    end
    if j==2
       Q(j,1)=u_w;
       Q(1,j)=0.5*u_w;
       Q(j,j+1)=0.5*u_w;
    elseif j>2
        Q(j,j-1)=0.5*u_w;
    end
end

% Constructing generalized eigenvalue problem in the form Dx=c*Fx
% Identity matrix
I = zeros(N+1);
for j=1:length(E)-4
   I(j,j) = 1; 
end
% Solving the generalized eigenvalue problem
D = (sqrt(-1)/(alpha*Re))*( E^4-2*(alpha^2)*I*(E^2)+(alpha^4)*I )+ 3*I - (alpha^2)*Q + Q*(E^2); % LHS
F = I*E^2 - (alpha^2)*I; % RHS
%BCs
for k = 1:length(D)
   D(N-2,k) = 1; % phi(1)=1
   D(N,k) = (k-1)^2; % phi'(1)=n^2
   D(N-1,k) = (-1)^(k-1); % phi(-1)=(-1)^(n-1)
   D(N+1,k) = ((-1)^(k-1))*(k-1)^2; % phi'(-1)= n^2 (-1)^(n-1)
end
% Solving for the eigenvalues and eigenvectors
format long
[V,W] = eig(D,F);
eigvals = diag(W);
eigvals = sort(eigvals) % Display eigenvalues
min(sort(diag(W))) % Display least stable eigenvalue

%% 

c = [];
Re = 3848;
alpha = 1.0;
u_w=0;
for m=9:9
    N = 2^m;
    gridpts = 2^(m-1); % # grid points
    h = 2/(gridpts-1);  
    U = zeros([gridpts 1]);
    for i = 1:gridpts-1
        U(i)= 1.5*(1-(-1 + h*(i-1))^2) + u_w*(-1 + h*(i-1));
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
        A(k,k+1) = 3 ; %phi_j
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













