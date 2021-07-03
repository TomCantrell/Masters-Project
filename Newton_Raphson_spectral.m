
% =================================================
% ---- Computation of critical Reynolds number ----
% =================================================

% Using the Newton-Raphson root-finding method to the critical Reynolds
% number, utilising the spectral implementation.

N=54; %N+1 vary the number of Chebyshev polynomials used in the approxn
Re=5500:50:7000;
alpha=1.02056;
Re_I = []; 
for j=1:length(Re)
    % Work out least stable eigenvalue
    c = Spectral(Re(j),alpha,N);
    if imag(c) > 0
        Re_I = [Re_I;Re(j-1)];
        Re_I = [Re_I;Re(j)];
        break
    end    
end
Re_I;

% %% fixed disturbance
correction = 1; % Initialise the correction
c_r = 2;
c=1;fprime=1;
alpha=1.02056;
while abs(correction) > 1e-10 % While the correction is greater than Re+-.5 continue to perform iterations until convergence
   Re_I(1) = Re_I(1) - imag(c)/fprime;
   fprime = (imag(Spectral(Re_I(1)+25,alpha,N))-imag(Spectral(Re_I(1)-25,alpha,N)))/(2*25);
   c = Spectral(Re_I(1),alpha,N)
   % Update estimate of root
   
   correction = imag(c);
   %correction = c_r - real(c);
   %c_r = real(c);
   Re_I(1);
end
Re_I(1)

%% 
% Performing the Newton-Raphson iterative algorithm to converge to the 
% critical Reynolds number
alpha=1.0202:0.00001:1.0208;
RE=[];
for k = 1:length(alpha)
    correction = 1; % Initialise the correction
    c_r = 2;
    while abs(correction) > 1e-11 % While the correction is greater than Re+-.5 continue to perform iterations until convergence
       fprime = (imag(Spectral(Re_I(1)+50,alpha(k),N))-imag(Spectral(Re_I(1)-50,alpha(k),N)))/(2*50);
       c = Spectral(Re_I(1),alpha(k),N)
       % Update estimate of root
       Re_I(1) = Re_I(1) - imag(c)/fprime;
       correction = imag(c);
       %correction = c_r - real(c);
       %c_r = real(c);
       Re_I(1);
    end
    RE = [RE;Re_I(1)];
    alpha(k);
end

min_Re = min(RE);
for i=1:length(RE)
    if RE(i)==min_Re
        t = i;
    end
end
min_Re
alpha(t)

plot(RE,alpha)
grid on







