%---------------- Linear 1D Eigenvalue Problem ------------------
n = 5;
h = 1/n;
A = zeros(n-1);
% create tridiagonal matrix
for i = 1:n-1
   for j = 1:n-1
       if i == j
          A(i,j) = -2;
          if i > n-1
              break
          else
              A(i+1,j) = 1;
              A(i,j+1) = 1;
          end
       end
    end
end







