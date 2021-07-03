
% ===============================================================
% Utilising a local iterative method to solve the generalzied
% eigenvalue problem more quickly.
% ===============================================================
% ( Boundary value problem )

% Guess f using f(z)=z(1-z) 
% 1D BVP 
lambdaG = 8.5;  % Initial guess for first eigenvalue (pi^2)  
for m=3:4
    N = 2^m + 1 
    h = 1/(N-1);
    fG=zeros([N 1]);
    % Using an arbitrary function that satisfies the boundary conditions
    for j=1:N
        fG(j) = (j-1)*h*(1-(j-1)*h) ;
    end
    o=1;
    l=[1];
    correction = 1;
    Correct=[];
    error=[abs(lambdaG-pi^2)];
    while abs(correction) > 10^-8
        N=2^m+1;
        A = zeros(N);
        B = zeros([N 1]);
        h = 1/(N-2);
        A(1,2) = -3/(2*h);  %phi'(0)=1 BC
        A(1,3) = 4/(2*h);
        A(1,4) = -1/(2*h);
        A(2,2) = 1;   % phi=0 BC
        B(1) = 1 - (-3*fG(1)+4*fG(2)-fG(3))/(2*h);
        B(2) = -fG(1);
        for j=3:N-1
            A(j,j) = lambdaG - 2/(h^2);
            A(j,j-1) = 1/(h^2);
            A(j,j+1) = 1/(h^2);
            A(j,1) = fG(j-1);
            B(j) = -(fG(j)-2*fG(j-1)+fG(j-2))/(h^2) - lambdaG*fG(j-1);
        end
        A(N,N) = 1; % phi=0 BC
        B;
        B(N) = -fG(N-1);
        
        % solve for the correction terms and update the guess
        solv = A\B;
        fGimprov = zeros([N 1]); 
        for k=2:N
             fGimprov(k-1) = solv(k);
        end
        fG = fG + fGimprov;
        correction = solv(1);
        Correct = [Correct; correction];
        lambdaG = lambdaG + solv(1);
        error=[error; abs(lambdaG-pi^2)];
        error
        o=o+1;
        l=[l;o];
    end
end


 

         