c = [];
Re = 5772.22;
alpha = 1.02056;
for m=9:9
    N = 2^m;
    gridpts = 2^(m-1); % # grid points
    h = 2/(gridpts-1);  
    U = zeros([gridpts 1]);
    for i = 1:gridpts-1
        U(i)= 1-(-1 + h*(i-1))^2;
    end
    A = zeros(N);
    B = zeros(N);
    b = (-1)^0.5 / (alpha*Re);
    % phi' = 0 BC
    A(1,2) = -3;
    A(1,4) = 4;
    A(1,6) = -1;
    % phi' = 0 BC
    A(2,2) = 1;
    for j = 2:N/2 - 1
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
    A(N-1,N-4) = 1;
    A(N-1,N-2) = -4;
    A(N-1,N) = 3;
    % phi = 0 BC
    A(N,N) = 1;
    [eigvec,v] = eig(A,B); % v generalized eigenvalues as entries on a diagonal matrix 
    vsorted = sort(diag(v));
    vsorted(1); % Least stable mode
    gridpts; % # grid points
    c = [c;abs(imag(vsorted(1)))];
end
v = diag(v);

% matrix of phi(z)'s
phi_eig_fns = zeros(N/2);
for j=1:N
    phi_temp = zeros([N/2 1]); % create temporary vector
    for k = 2:2:N
        phi_temp(k/2) = eigvec(k,j);
    end
    phi_eig_fns(:,j) = phi_temp;
end

imagc = [];
realc = [];
vec1 = []; vec2 = []; vec3 = [];
for j=1:N
    if imag(v(j)) > -1
        imagc = [imagc; imag(v(j))];  
        realc = [realc; real(v(j))];
        if real(v(j)) < 0.65
            vec1 = [vec1;j];
        end
        if imag(v(j)) < -0.4
            vec2 = [vec2; j];
        end
        if real(v(j)) > 0.7
            vec3 = [vec3;j];
        end     
    end
end
z = linspace(-1,1,N/2);
% Creating plots of phis
% redefining vecs
for k=1:length(vec1)
   vec1(k) = v(vec1(k));
end
minv1 = min(vec1);
figure(4)
for i=1:length(v)
   if v(i)==minv1
      i
      v(i)
      plot(real(phi_eig_fns(:,i)),z,"color","k")
      hold on
      plot(imag(phi_eig_fns(:,i)),z,"-.","color","b")
      grid on
      xlabel("\phi")
      ylabel("z")
   end
end
legend("Re(\phi)", "Im(\phi)", "Location", "best")
figure(2)
for k=1:length(vec2)
   vec2(k) = v(vec2(k));
end
minv2 = min(vec2);
for i=1:length(v)
   if v(i)==minv2
      i
      v(i)
      plot(real(phi_eig_fns(:,i)),z,"color", "k")
      hold on
      plot(imag(phi_eig_fns(:,i)),z, "-.","color","b")
      grid on
      xlabel("\phi")
      ylabel("z")
   end
end
legend("Re(\phi)", "Im(\phi)", "Location", "best")
figure(3)
for k=1:length(vec3)
   vec3(k) = v(vec3(k));
end
minv3 = min(vec3);
for i=1:length(v)
   if v(i)==minv3
      i
      v(i)
      plot(real(phi_eig_fns(:,i)),z,"color","k")
      hold on
      plot(imag(phi_eig_fns(:,i)),z, "-.","color","b")
      grid on
      xlabel("\phi")
      ylabel("z")
   end
end
legend("Re(\phi)", "Im(\phi)", "Location", "best")


% Finding critical phi
for i=1:N
    if v(i) == vsorted(1)
        phi_crit =[]; % phi(z) when c = c_crit
        for j=2:2:N
            phi_crit = [phi_crit; eigvec(j,i)];
        end
    end
end

% What do phi(z) fns look like?

% phi_crit
% figure(4)
% z = linspace(-1,1,length(phi_crit));
% plot(z,real(phi_crit))
% hold on
% plot(z,imag(phi_crit),"--")


phi_bound = [];
for k=1:length(v)-2 % -2 as final two eigenfns are just one at each end
    %eigvec(:,j)
    phi = [];
    for j=2:2:N
        phi = [phi; eigvec(j,k)];
        if eigvec(j,k) == 1 % trying to find where eigenfn is just 1
            phi_bound = [phi_bound;k];
        end
    end
    %plot(real(phi))
    %hold on
    %grid on
end

% plot of 5 least stable modes
% figure(5) 
% z = linspace(-1,1,N/2);
% for j=1:5
%     for i=1:N
%         if v(i) == vsorted(j)
%             phi_least_stable = [];
%             for j=2:2:N
%                 phi_least_stable = [phi_least_stable; eigvec(j,i)];
%             end
%             plot(z,real(phi_least_stable))
%             hold on
%         end
%     end
% end
% grid on
% %title("5 Least stable modes")
% ylabel("\phi(z)")
% xlabel("z")

phiz = phi_crit;
M = length(phiz);
% phi'(z)=0
phi_z = zeros([M 1]);
phi_z(1) = (-3*phiz(1)+4*phiz(2)-phiz(3))/(2*h);
phi_z(N/2) = (3*phiz(M)-4*phiz(M-1)+phiz(M-2))/(2*h);
for k=2:M - 1
    phi_z(k) = (phiz(k+1)-phiz(k-1))/(2*h);
    k;
end
c = vsorted(1);
c= real(vsorted(1)); % setting c_i = 0
x = linspace(0,2*pi/alpha,M);
z = linspace(-1,1,M);
U_hat = zeros(M);
W_hat = zeros(M);
t=0;
for k=1:M % varying x only for one cycle
    U_hat(:,k) = phi_z * exp( sqrt(-1)*alpha*(x(k)-c*t));% phi'(z)
    W_hat(:,k) = phiz*(-sqrt(-1))*alpha*exp( sqrt(-1)*alpha*(x(k)-c*t));% -i*alpha*phi(z)
end
% Create new array so that quiver plot is more readable
approx = 64;
U_hat_approx = zeros(M/approx); W_hat_approx = zeros(M/approx);
for l=1:M/approx
    for k=1:M/approx
        % fixed column, go down the column in 16's
        U_hat_approx(l,k) = U_hat(l*approx,k*approx);
        W_hat_approx(l,k) = W_hat(l*approx,k*approx);
    end 
end

figure(7)
% Redefine x, z for plot
x = linspace(0,2*pi/alpha,M/approx);
z = linspace(-1,1,M/approx);
quiver(x,z,real(U_hat_approx),real(W_hat_approx),"k");
grid on
xlabel("x")
ylabel("z")

% contour(x,z,real(U_hat_approx),real(W_hat_approx));
% grid on 
% hold off
% xlabel("x")
% ylabel("z")

% 
% for j=1:M
%         w = mode(x(j),t,alpha,c);
%         uhat = [uhat; phi_z(k)*w]; % phi'(z)
%         what = [what; -sqrt(-1)*alpha*phiz(k)*w]; % -i*alpha*phi(z)
%         t=t+30;
%     end



% figure(3)
% for i=1:M
%     plot(x,real(Uhat(:,i)),z,real(What(:,i))); % real(What(:,i))
% end



