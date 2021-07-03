
% Function that computes the least stable eigenvalue for given wavenumber
% Reynolds number and wall speed in Poiseuille-Couette flow by using the 
% spectral-tau method already developed but modified slightly for the
% hybrid flow
% ----------------------------------------------------------------
function c = spectral_PPCF(Re,alpha,u_w,N)
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
    [~,W] = eig(D,F);
    eigvals = diag(W);
    eigvals = sort(eigvals); % Display eigenvalues
    
    % Sanity check
    eig_crit = min(eigvals);
    for k=1:length(eigvals)-6
        if imag(eigvals(k)) > imag(eig_crit)
            eig_crit = eigvals(k);
        end
    end
    c = eig_crit;
end













