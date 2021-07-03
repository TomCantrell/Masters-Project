
% Spectral implementation of harmonic problem
format long
%%
N=8;
M = N+2;
A = zeros(M);
B = zeros(M);
% Chebyshev points, z ranges from 0 to 1
k = 0:1/M:1; z = zeros([length(k) 1]);
for i=1:length(k)
    z(i) = cos((i-1)*pi/(2*M));
end
% Here we set up the generalized eigenvalue problem
%BCs f(0)=0
% Fill matrices

% BC f(0)=0 
for k=1:M
    A(1,k)= ((1+(-1)^(k-1))/2) * (-1)^((k-1)/2);
end
for j=2:M-1
    n = j-2; %n runs from 0 to N+2
    for p=n+2:M-1
       if mod(p-n,2)==0
          A(j,p+1) = p*(p^2-n^2);
       end
    end
    if n==0
        B(j,1)=-2;
    else
        B(j,j-1)=-1;
    end
end
% BC f(1)=0
for k=1:M
   A(M,k) = 1;
end
A;
% Solve the eigenvlaue problem
[V,W] = eig(A,B);
v = sort(diag(W))
%%
format long

N=8;
A=zeros(N);
B=zeros(N);
C=zeros(N);

% Constructing C matrix to be inverted and then right multiplied by matrix
% A to find a matrix E that can be used recursively to find derivatives
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
E = C\A


% Constructing E matrix directly 
E=zeros(N+1);
for n=1:N
   for m=n:2:N % n+1 but n here bc indexing
      if n-1 == 0
          c=2;
      else
          c=1;
      end
      E(n,m+1) = (2/c) * m; 
   end   
end
E
A=4*E*E; % 2nd derivative 
% Add BCs
for j=1:N+1
   A(N,j) = (-1)^(j-1);
   A(N+1,j)=1;
end
% C matrix Ax = lambda Cx
C = zeros(N+1);
for j=1:N-1
   C(j,j)=-1; 
end
C;
% Working out eigenvalues
[V,W] = eig(A,C);
sort(diag(W))







