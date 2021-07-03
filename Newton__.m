
% ==================================================================
% --  Local iterative method used in obtaining the neutral curve  --
% ==================================================================

%

function [PsiG] = Newton__(PsiG,Re,alpha)
    correction=[];
    cG = PsiG(1);
    N = length(PsiG)-1;
    h = 2/((N/2)-1);
    % U 
    U = zeros([N/2 1]);
    for i = 1:N/2-1
        U(i)= 1-(-1 + h*(i-1))^2;
    end
    %=======================================
    c_imag=1;
    x=0; 
    while c_imag > 1e-8
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
            b(j) = -((U(i)-cG-sqrt(-1)*(2/h^2 +alpha^2)/(alpha*Re))*PsiG(j) + 2*PsiG(j+1) + sqrt(-1)*(PsiG(j+2)+PsiG(j-2))/(alpha*Re*h^2));
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
        c = PsiG(1);
    end 
end














