
% =================================================
% ---- Computation of critical Reynolds number ----
% =================================================

% Newton-Raphson implementation to obtain an approximation to the critical
% Reynolds number in plane Poiseuille flow, by using finite-difference
% methods.

% Our initial guess required
Re = 5750;
alpha = 1.02056;
m=11;
n = 2^m;
t = 2^(m-1); % # grid points
h = 2/(t-1);  
U = zeros([t 1]);
for i = 1:t-1
    U(i)= 1-(-1 + h*(i-1))^2;
end
A = zeros(n);
B = zeros(n);
b_ = sqrt(-1)/(alpha*Re);
% phi' = 0 BC
A(1,2) = -3/(2*h);
A(1,4) = 4/(2*h);
A(1,6) = -1/(2*h);
% phi = 0 BC
A(2,2) = 1;
for j = 2:n/2 - 1
    k = 2*j-1;
    B(k,k) = 1;
    A(k,k-2) = b_/h^2 ; %psi_j-1
    A(k,k) = U(j) - (2/h^2 + alpha^2)*b_; %psi_j
    A(k,k+1) = 2 ; %phi_j
    A(k,k+2) = b_/h^2; %psi_j+1
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
%% ------------------------------------------------------------------------
%------------------------------------------------------------------------
% interpolation of eigenvector to approximate to a larger grid
% interpolate until size of Psiginter exceeds 
q=10;
Psiginter = 1;
while length(Psiginter) < 2^q
    Psiginter = zeros([2*(length(Psig)) 1]);
    i=1;
    for k=1:length(Psig)/2   
        if k >=length(Psig)/4 + 1 
            l=2*(k-1)+1;
            Psiginter(2*l+1)=  Psig(i); % Psi eqn
            Psiginter(2*l+2)=  Psig(i+1); % Phi eqn
            i=i+2;
        else
            l=2*(k-1)+1;
            Psiginter(2*l-1)=  Psig(i); % Psi eqn
            Psiginter(2*l)=  Psig(i+1); % Phi eqn
            i=i+2;
        end
    end
    p=1;
    for j=1:length(Psig)/2-1
        if j==length(Psig)/4
            r=4*j-1;
            Psiginter(r)= (Psig(p)+Psig(p+2))/2; %
            Psiginter(r+1)= (Psig(p+1)+Psig(p+3))/2; %
            Psiginter(r+2)= (Psig(p)+Psig(p+2))/2;
            Psiginter(r+3)= (Psig(p+1)+Psig(p+3))/2;
            p=p+2;
        elseif j>length(Psig)/4
            r = 4*j+1;
            Psiginter(r)= (Psig(p)+Psig(p+2))/2; %
            Psiginter(r+1)= (Psig(p+1)+Psig(p+3))/2; %            
            p=p+2;
        else 
            r = 4*j-1;
            Psiginter(r)= (Psig(p)+Psig(p+2))/2; %
            Psiginter(r+1)= (Psig(p+1)+Psig(p+3))/2; %
            p=p+2;
        end
    end
    Psig = Psiginter;
    n = length(Psiginter);
end
N = length(Psiginter);
PsiG = zeros([N+1 1]);
PsiG(1) = cG;
for j=2:N+1
    PsiG(j) = Psiginter(j-1);
end
% ---------------------------------------------------------------------
%% ------------------------------------------------------------------------
% Now, implementing the Newton-Raphson root-finding method

% Finding where the sign changes i.e. an interval of Reynolds numbers where
% the root will lie
Re=5750:50:6050;
Re_I = []; 
for j=1:length(Re)
    % Work out least stable eigenvalue
    v = Newton__(PsiG,Re(j),alpha);
    v(1)
    if imag(v(1)) > 0
        Re_I = [Re_I;Re(j-1)];
        Re_I = [Re_I;Re(j)];
        break
    end    
    PsiG=v;
end
Re_I
correction=10;
%%

tic
while abs(correction)>0.1
   % Performing Newton-Raphson iteration
   u = Newton__(PsiG,Re_I(1)+25,alpha);
   w = Newton__(PsiG,Re_I(1)-25,alpha);
   fprime = (imag(u(1))-imag(w(1)))/(2*25);
   c_i = imag(v(1))
   % Update estimate of root
   Re_I(1) = Re_I(1) - c_i/fprime
   correction = c_i/fprime;  
   v = Newton__(PsiG,Re_I(1),alpha);
   PsiG = v;
   imag(PsiG(1)); % this should --> 0
end
toc
Re_I(1)



%% 
% want to extrapolate to infinity
Re_crit = [5821.89,5784.81,5775.4,5773.02];
dcrmnt = [];
for k=1:length(Re_crit)-2
    dcrmnt = [dcrmnt;(Re_crit(k)-Re_crit(k+1))/(Re_crit(k+1)-Re_crit(k+2))];
end
Re_G = [5775.4,5773.02];
e = Re_crit(k+1)-Re_crit(k+2);
Re_g = Re_G(2);
while e>1e-5
    Re_g = Re_g - e/4;
    %Re_crit = [Re_crit;Re_G];
    Re_G(1) = Re_G(2);
    Re_G(2) = Re_g;
    Re_g
    e = Re_G(1)-Re_G(2);
end












