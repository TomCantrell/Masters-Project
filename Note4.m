% % 

 %%%   SAVED 22/12/20 10:55
% % =====================  Newton's Iteration ============================
% % A function that takes guess and then iterates
% function [Psig] = Newton(Re,alpha,PsiG)
% ===================== Inital Guess ====================
c_r = [];
Re = 5772.22;
alpha = 1.02056;
m=7;
n = 2^m;
t = 2^(m-1); % # grid points
h = 2/(t-1);  
U = zeros([t 1]);
for i = 1:t-1
    U(i)= 1-(-1 + h*(i-1))^2;
end
A = zeros(n);
B = zeros(n);
b = sqrt(-1)/(alpha*Re);
% phi' = 0 BC
A(1,2) = -3/(2*h);
A(1,4) = 4/(2*h);
A(1,6) = -1/(2*h);
% phi = 0 BC
A(2,2) = 1;
for j = 2:n/2 - 1
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
A(n-1,n-4) = 1/(2*h);
A(n-1,n-2) = -4/(2*h);
A(n-1,n) = 3/(2*h);
% phi = 0 BC
A(n,n) = 1;
[V,D] = eig(A,B);
eigsort = sort(diag(D));
v = eigsort;
cG = v(1);  
for i=1:n
    if D(i,i) == cG
        k=i;
    end
end
Psig = V(:,k:k); % Obtains eigenvector from corresponding eigenvalue
%------------------------------------------------------------------------
% interpolation of eigenvector to approximate to a larger grid
% interpolate until size of Psiginter exceeds 
o=1;
Psiginter = 1;
while length(Psiginter) < 1500
    Psiginter = zeros([2*(length(Psig)/2 -1)+n 1]);
    i=1;
    for k=1:4:length(Psiginter)    
        Psiginter(k)=  Psig(i); % Psi eqn
        Psiginter(k+1)=  Psig(i+1); % Phi eqn
        i=i+2;
    end
    p=1;
    for j=3:4:length(Psiginter)
        Psiginter(j)= (Psig(p)+Psig(p+2))/2; %
        Psiginter(j+1)= (Psig(p+1)+Psig(p+3))/2; %
        p=p+2;
    end

    Psig = Psiginter;
    n = length(Psiginter);
    o=o+1;
end
N = length(Psiginter);
PsiG = zeros([N+1 1]);
PsiG(1) = cG;
for j=2:N+1
    PsiG(j) = Psiginter(j-1);
end

correction=[];
cG = PsiG(1); 
N = length(PsiG)-1;
h = 2/((N/2)-2);
% U 
U = zeros([N/2 1]);
for i = 1:N/2-1
    U(i)= 1-(-1 + h*(i-1))^2;
end

c_imag=1;
x=0; 
while c_imag > 10^-8
    A = zeros(N+1);
    b = zeros([N+1 1]);
    %----------BCs--------------
    % psi = 1 BC
    A(1,2) = 1;
    % phi' = 0 BC
    A(2,3) = -3/(2*h);
    A(2,5) = 4/(2*h);
    A(2,7) = -1/(2*h);
    % phi = 0 BC
    A(3,3) = 1;
    % First 3 rows of B vector
    b(1) = 1-PsiG(2); % psi=1 BC
    b(2) = -(-3*PsiG(3)+4*PsiG(5)-PsiG(7))/(2*h); % phi'=0 BC
    b(3) = -PsiG(3); % phi=0 BC 
    %--------------------------
    for i = 2:N/2 - 1
        j = 2*i;
        cG = PsiG(1);
        % B vector
        b(j) = -( (U(i)-cG-sqrt(-1)*(2/h^2 +alpha^2)/(alpha*Re))*PsiG(j) + 2*PsiG(j+1) + sqrt(-1)*(PsiG(j+2)+PsiG(j-2))/(alpha*Re*h^2));
        b(j+1)= -( PsiG(j) - (PsiG(j+3)+PsiG(j-1))/h^2 + (alpha^2 + 2/(h^2))*PsiG(j+1));
        % A matrix
        % Fill psi eqn
        A(j,j) = U(i)- cG - sqrt(-1)/(alpha*Re) * (alpha^2 + 2/h^2); % psi_j
        A(j,1) = -PsiG(j); % c-tilde
        A(j,j+1)= 2 ; % phi_j
        A(j,j-2) = sqrt(-1)/(alpha*Re*h^2) ; %psi_j-1
        A(j,j+2) = sqrt(-1)/(alpha*Re*h^2); % psi_j+1
        % Fill phi eqn
        A(j+1,j+1) = 2/(h^2) + alpha^2 ; % phi_j
        A(j+1,j) = 1 ; % psi_j
        A(j+1,j-1) = -1/(h^2); % phi_j-1
        A(j+1,j+3) = -1/(h^2); % phi_j+1
    end
    %--------BCs------------
    % Last two rows of B vector 
    b(N)= -(3*PsiG(N+1)-4*PsiG(N-1)+PsiG(N-3))/(2*h);
    b(N+1)= -PsiG(N+1);
    % phi' = 0 BC
    A(N,N-3) = 1/(2*h);
    A(N,N-1) = -4/(2*h);
    A(N,N+1) = 3/(2*h);
    % phi = 0 BC
    A(N+1,N+1) = 1;
    %----------------------
    % Solve for correction and update guess
    solv = A\b;
    PsiG = PsiG + solv;

    c_imag = abs(imag(solv(1)));

    x=x+1;
    correction = [correction; c_imag]; % correction c~
    c = PsiG(1)
end
Psig = PsiG;
%end

