
% =====================================================================
% ----- Spectral method implementation - Orr-Sommerfeld equation ------
% =====================================================================


%% Formulate the generalized eigenvalue problem 
Re = 5772.22; % Reynolds number
alpha = 1.02056; % Wavenumber
N=119; % N+1 Chebyshev polynomials
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
        Q(j,j+2) = -0.25;
    end
    % Diagonal entries
    if j-1 == 1
        Q(j,j)=0.25;
    else
        Q(j,j)=0.5;
    end        
    % Other off-diagonal entries
    if j>2
        if j==3
            Q(j,j-2) = -0.5;
        else
            Q(j,j-2)=-0.25;
        end
    end
end
% Constructing generalized eigenvalue problem in the form Dx=cFx
% Identity matrix
I = zeros(N+1);
for j=1:length(E)-4
   I(j,j) = 1; 
end
% Solving the generalized eigenvalue problem
D = (sqrt(-1)/(alpha*Re))*( E^4-2*(alpha^2)*I*(E^2)+(alpha^4)*I )+ 2*I - (alpha^2)*Q + Q*(E^2); % LHS
F = I*E^2 - (alpha^2)*I; % RHS
%Boundary conditions
for k = 1:length(D)
   D(N-2,k) = 1; % phi(1)=1
   D(N,k) = (k-1)^2; % phi'(1)=n^2
   D(N-1,k) = (-1)^(k-1); % phi(-1)=(-1)^(n-1)
   D(N+1,k) = ((-1)^(k-1))*(k-1)^2; % phi'(-1)= n^2 (-1)^(n-1)
end
% ----------   Solving for the eigenvalues and eigenvectors  --------------
% Columns of V are the eigenvectors, with W a diagonal matrix with
% eigenvalues as the entries
format long
[V,W] = eig(D,F); % 
eigvals = diag(W);
eigvals = sort(eigvals) % Display eigenvalues

% Sanity check
eig_crit = min(eigvals);
for k=1:length(eigvals)-6 %-6 as we need to remove the spurious eigenvalues
    if imag(eigvals(k)) > imag(eig_crit)
        eig_crit = eigvals(k);
    end
end
eig_crit % display least stable eigenvalue

for j=1:length(eigvals)
   if W(j,j) == eig_crit
       eigvec_crit = V(1:length(eigvals),j);
   end
end
% Approximation of coefficients, eigvec_crit has coefficients equal to zero
% for the odd order Chebyshev polynomials 

% Want to find positions of 10 least stable eigenvalues 
f = zeros([10 1]);
s = zeros([length(f) 1]); ctr=0; eigvals_im = flip(sort(imag(eigvals)));
for k=5:length(f)+4
    f(k-4) = eigvals_im(k);
end
%Positions of least stable eigenvalues
for i=1:length(f)
    ctr=ctr+1;
    for k=1:N+1
        if imag(W(k,k)) == f(ctr)
            imag(W(k,k));
            eigvals_im(ctr);
            s(ctr) = k;
        end
    end
end

%-------   Plotting the data obtained from solving the OS problem  --------
% Plot eigenvalue spectrum
close(figure(1))
figure(1)
plot(eigvals(1:length(eigvals)),"o","color","k")
ylim([-1 0.1]);
xlim([0.2 1]);
grid on
hold on
plot(real(eig_crit),imag(eig_crit),"*","color","r","LineWidth",2)
xlabel("c_r")
ylabel("c_i")

%% Plot eigenfunctions and vector field
figure(2)
% Generate an approximation to the eigenfunction using the coefficients
% solved for in the generalized eigenvalue problem
z=-1:0.005:1; eig_fn=zeros([1 length(z)]); eig_fn_1=zeros([1 length(z)]);  % initialise vectors

% Other eigenvalue/vector
p = 69;V(1:N+1,p);% % 94 N=120
coeff = eigvec_crit; % Vector of coefficients corresponding to least stable eigenvalue
coeff_1 = E*coeff; % Coefficients of the 1st derivative
for k=1:N+1
    eig_fn = eig_fn+coeff(k)*chebyshevT(k-1,z); % Linear combination of Chebyshev polynomials  
    eig_fn_1 = eig_fn_1 + coeff_1(k)*chebyshevT(k-1,z);
end

plot(real(eig_fn),z,"LineWidth",0.9,"color","k")
hold on
plot(imag(eig_fn),z,"--","LineWidth",0.9,"color","b")
xlabel("\phi")
ylabel("z")
legend("Re(\phi)","Im(\phi)","Location","Best")
grid on

% Vector field plot + contours
c = eig_crit;
%c = W(p,p)
x = linspace(0,6*pi/alpha,length(z));
u = zeros(length(eig_fn)); 
w = zeros(length(eig_fn));

t=0;
for j=1:length(eig_fn)
   u(1:length(eig_fn),j) = eig_fn_1*exp(sqrt(-1)*alpha*(x(j)-c*t));
   w(1:length(eig_fn),j) = -sqrt(-1)*alpha*eig_fn*exp(sqrt(-1)*alpha*(x(j)-c*t));
end

l=0;
approx = 6; % this controls the number of arrows plotted i.e resolution of vector field
x_approx = []; z_approx = [];
for k=1:approx:length(z)
    l = l+1;
    x_approx(l) = x(k);
    z_approx(l) = z(k);
end

m=0;
l = length(x_approx);
u_approx = zeros(l);
w_approx = zeros(l);
for k=1:approx:length(z)
    m = m+1;
    i=0;
    for j=1:approx:length(z)
        i=1+i;
        u_approx(i,m) = u(j,k);
        w_approx(i,m) = w(j,k);
    end
end
% Plotting vector field
% Contour of magnitude of vectors
X = zeros(length(w));
for k=1:length(X)
    for j=1:length(X)
       X(j,k) = sqrt((u(j,k)^2+w(j,k)^2)); 
    end
end

close(figure(3))
figure(3)
contourf(x,z,real(X),"LineWidth",0.25)
hold on
%q=quiver(x_approx,z_approx,real(u_approx),real(w_approx),"k","LineWidth",0.9);
% q.ShowArrowHead = 'off';
% q.Marker = '.';
xlabel("x")
ylabel("z")
hold on
%yline(1,"color","k","LineWidth",2)
hold on
%yline(-1,"color","k","LineWidth",2)
ylim([-1 1])
xlim([0 4*pi])
grid on

close(figure(4))
figure(4)
contourf(x,z,real(u),12,"LineWidth",0.8)
xlabel("x")
ylabel("z")
grid on









