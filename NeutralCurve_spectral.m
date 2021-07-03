
% ===============================================================
% -----  Linear neutral curve obtained by spectral methods  ----- 
% ===============================================================

%%
% 500secs
Re_upper=[];
Re_lower=[];
a_upper=[];
a_lower=[];
c_r_upper=[];
c_r_lower=[];
%alpha=1.09; % Wavenumber
a=0.50:0.005:1.065; % Wavenumber alpha
for p=1:length(a)
    N=74; % N+1 Chebyshev polynomials
    % Generate the approximation to the function f
    Re=5600:300:60000;
    c_i=length(Re);
    for k=1:length(Re)
        c_i(k)= imag(Spectral(Re(k),a(p),N));
    end
    % ===========   Newton-Raphson root-finding method   ============
    % This loop locates a sign change from the vector generated above on a
    % coarse grid. These values are then used as the guesses for the
    % root-finding method. 
    c_i_IG = []; % Vector of initial guesses (IG)
    Re_IG = []; 
    for j=2:length(c_i)
        if c_i(j)>0 && c_i(j-1)<0
            c_i_IG = [c_i_IG,c_i(j-1)];
            Re_IG = [Re_IG,Re(j-1)];
        end
        if (c_i(j)<0) && (c_i(j-1)>0)
            c_i_IG = [c_i_IG,c_i(j-1)];
            Re_IG = [Re_IG,Re(j-1)];
        end
        if length(Re_IG)==2
            break
        end
    end
    if isempty(Re_IG)
       % Move on
    else
        for k=1:length(Re_IG)
            correction = 10; % Initialise the correction
            while abs(correction) > 1 % While the correction is greater than Re+-.5 continue to perform iterations until convergence
               fprime = (imag(Spectral(Re_IG(k)+25,a(p),N))-imag(Spectral(Re_IG(k)-25,a(p),N)))/(2*25);
               c_i = imag(Spectral(Re_IG(k),a(p),N));
               % Update estimate of root
               Re_IG(k) = Re_IG(k) - c_i/fprime;
               correction = c_i/fprime;
            end
        end
    end
    if length(Re_IG)>1
        Re_lower = [Re_lower;Re_IG(1)];
        Re_upper = [Re_upper;Re_IG(2)];
        a_lower = [a_lower;a(p)];
        a_upper = [a_upper;a(p)];
        c_r_upper = [c_r_upper;real(Spectral(Re_IG(2),a(p),N))];
        c_r_lower = [c_r_lower;real(Spectral(Re_IG(1),a(p),N))];
    elseif length(Re_IG)==1
        Re_lower = [Re_lower;Re_IG(1)];
        a_lower = [a_lower;a(p)];
        c_r_lower = [c_r_lower;real(Spectral(Re_IG(1),a(p),N))];        
    end
    Re_IG
end
%% Finer part of the upper branch
N=74;
Re_u = 5800:100:20000; % Reynolds numbers corresponding to the upper part of the neutral curve
alpha_u=[];
c_r_u = [];
for p=1:length(Re_u)
    alpha=0.85:0.0025:1.15;
    c_i_upper=zeros([length(alpha) 1]);
    for j=1:length(alpha)
        c_i_upper(j) = imag(Spectral(Re_u(p),alpha(length(alpha)-j+1),N));
        %alpha(length(alpha)-j+1)
    end
    for k=2:length(c_i_upper)
       if c_i_upper(k)>0 && c_i_upper(k-1)<0 
          %Re_u = [Re_u;Re];
          alpha_g = alpha(length(alpha)-k+1); % guess for alpha
          break
       end
    end
    alpha_g;
    % =======  Newton-Raphson implementation for the upper branch  =======
    correction=10;
    while correction>10^(-6)
       % Perform iterations of Newton-Raphson root-finding algorithm until the threshold in the while loop is reached
       fprime = (imag(Spectral(Re_u(p),alpha_g+0.05,N))-imag(Spectral(Re_u(p),alpha_g-0.05,N)))/(2*0.05);
       c_i = imag(Spectral(Re_u(p),alpha_g,N));
       % Update estimate of root
       alpha_g = alpha_g - c_i/fprime;
       correction = c_i/fprime;
    end
    c_r_u = [c_r_u;real(Spectral(Re_u(p),alpha_g,N))];
    alpha_g
    Re_u(p);
    alpha_u = [alpha_u;alpha_g];
end
alpha_u;
%% Plotting the neutral curve with the curve for wavespeed
figure(1)
plot(Re_lower,a_lower,"color","k","LineWidth", 1.5)
hold on
plot(Re_upper,a_upper,"color","k","LineWidth",1.5)
hold on
plot(Re_u,alpha_u,"color","k","LineWidth",1.5)
xlabel("Re")
ylabel("\alpha")
grid on
xlim([0 50000])
ylim([0 1.2])

figure(2)
plot(Re_upper,c_r_upper,"color","k","LineWidth",1.5)
hold on
plot(Re_u,c_r_u,"color","k","LineWidth",1.5)
hold on
plot(Re_lower,c_r_lower,"color","k","LineWidth",1.5)
grid on
xlabel("Re")
ylabel("c_r")
%legend("Im(c)=0")
ylim([0 0.3])
xlim([0 50000])


















































































%%

% % Brute force approach of trying every point
% 
% %% Formulate the generalized eigenvalue problem 
% 
% % Working out E matrix 
% N=98; % 98 is a good number
% % Define vectors for alpha and Re to loop over
% % Resolution (Re, alpha) = (5, 0.05)
% 
% % Create matrix with imaginary part of least stable eigenvalue for a given
% % set of (Re, alpha) values 
% alpha = [];
% for a=0.8:0.01:1.15
%    alpha = [alpha;a]; 
% end
% Re = [];
% for r=5500:20:8000
%    Re = [Re;r]; 
% end
% c_i = zeros([length(alpha) length(Re)]);
% c_r = zeros([length(alpha) length(Re)]);
% A = zeros(N);
% B = zeros(N);
% C = zeros(N);
% for i=1:N
%     if i==1
%         C(i,i)=2;
%     else 
%         C(i,i)=1;
%     end
%     if i<=N-2
%         C(i,i+2)=-1;  
%     end
%     A(i,i+1)=2*i;
% end
% A; 
% C; 
% E = C\A;
% E(N+1,N)=0; % Add another row to make the matrix square
% % Constructing the matrix which is used to determine the matrix formulation
% % of the cross terms in the Orr-Sommerfeld equation U*phi, U*phi"
% Q = zeros(N+1); % has same dimensions as E
% for j=1:N-3
%     if j<N
%         Q(j,j+2) = -0.25;
%     end
%     % Diagonal entries
%     if j-1 == 1
%         Q(j,j)=0.25;
%     else
%         Q(j,j)=0.5;
%     end        
%     % Other off-diagonal entries
%     if j>2
%         if j==3
%             Q(j,j-2) = -0.5;
%         else
%             Q(j,j-2)=-0.25;
%         end
%     end
% end
% % Constructing generalized eigenvalue problem in the form Dx=c*Fx
% % Identity matrix
% I = zeros(N+1);
% for j=1:length(E)-4
%    I(j,j) = 1; 
% end
% % ---  Solving the generalized eigenvalue problem  ---
% 
% for i=1:length(Re)
%     Re(i)
%     for j=1:length(alpha)
%         
%         D = (sqrt(-1)/(alpha(j)*Re(i)))*( E^4-2*(alpha(j)^2)*I*(E^2)+(alpha(j)^4)*I )+ 2*I - (alpha(j)^2)*Q + Q*(E^2);
%         F = I*E^2 - (alpha(j)^2)*I;
%         %BCs
%         for k = 1:length(D)
%            D(N-2,k) = 1; % phi(1)=1
%            D(N,k) = (k-1)^2; % phi'(1)=n^2
%            D(N-1,k) = (-1)^(k-1); % phi(-1)=(-1)^(n-1)
%            D(N+1,k) = ((-1)^(k-1))*(k-1)^2; % phi'(-1)= n^2 (-1)^(n-1)
%         end
%         % Solving for the eigenvalues and eigenvectors
%         [V,W] = eig(D,F);
%         eigvals = diag(W);
%         
%         
%         
%         c_r(length(alpha)-j+1,i) = real(min(eigvals));
%         c_i(length(alpha)-j+1,i) = imag(min(eigvals)); 
%     end
% end