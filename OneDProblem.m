%---------------- Linear 1D Eigenvalue Problem ------------------
n=100;
k = 1; % First k Eigenvalues 
eig1 = [];
error1 = [];
H = [];
for N=1:1:n
    h = 1/N;
    A = zeros(N);
    B = zeros(N);
    z = 0:1/N-1:1; % creating tridiagonal matrix
    A(1,1) = -1/h^2;
    A(N,N) = -1/h^2;
    for i = 2:N-1
        for j = 2:N-1
            if i==j
                B(i,j) = 1;
                A(i,j)= 2/h^2;
                if i+1 < N
                    A(i,j+1) = -1/h^2;
                    A(i+1,j) = -1/h^2;
                end
            end
        end
    end
    v = eig(A,B); % Eigenvalues not in order
    v = sort(v);
    
    Eig_exact = zeros([k 1]);
    m = zeros([k N]);
    for l = 1:1:k
        Eig_exact(l) = (l^2)*(pi^2);
    end
    if length(v) >= 3
        eig1 = [eig1;v(3)];
        error1 = [error1 ; abs(Eig_exact - eig1(N-2))];
        H = [H; h];
    end
end
H = log(H);
error1 = log(error1);

plot(H, error1)
r = (error1(1)-error1(N-2)) / (H(1)-H(N-2)) ;



