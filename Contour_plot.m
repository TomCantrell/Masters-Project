% DOESNT REALLY WORK
%  ================ Contour plot of the Neutral curve  ====================

% ===================== Inital Guess ====================
c_r = [];
Re = 5780;
alpha = 1.02;
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
while length(Psiginter) < 1000
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
% =========================================================================

Re = 5500:200:20500;
alpha = 0.6:0.05:1.15;
c = zeros([length(alpha) length(Re)]);
col=1;
for Re = 5500:200:20500
    c_tmp = [];
    for alpha = 0.6:0.05:1.15
        v = Newton(Re,alpha,PsiG);
        c_tmp = [c_tmp; v(1)];
        PsiG=v;
    end
    c(:,col) = c_tmp;
    col=col+1;
    Re
end
contour(Re,alpha,imag(c))

%  ~ 5 mins












