% Needs a redo

% ======================================================================
% ---------------------   Plot of convergence   ------------------------
% ======================================================================
Re = 5772.22;
alpha = 1.02056;
c = [];gridpoints=[];

for m=5:10
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
    c = [c;abs(imag(cG))];
    gridpoints = [gridpoints;t];
end

%% ------------------------------------------------------------------------
%------------------------------------------------------------------------
% interpolation of eigenvector to approximate to a larger grid
% interpolate until size of Psiginter exceeds 
for q=11:13
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
    
    %
    PsiG = Newton__(PsiG,Re,alpha);cG=PsiG(1);
    c=[c;abs(imag(PsiG(1)))];
    gridpoints = [gridpoints;2^(q-1)];
end
gridpoints=[gridpoints;2^13];
c = [c;0.00000132];
%% Plotting the data
plot(log(gridpoints),log(c))



