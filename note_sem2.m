
%
%%
N=80;
a=1.02056;
Re=5600:300:50000;
c_i=length(Re);
for k=1:length(Re)
    c_i(k)= imag(Spectral(Re(k),a,N));
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
           fprime = (imag(Spectral(Re_IG(k)+25,a,N))-imag(Spectral(Re_IG(k)-25,a,N)))/(2*25);
           c_i = imag(Spectral(Re_IG(k),a,N));
           % Update estimate of root
           Re_IG(k) = Re_IG(k) - c_i/fprime
           correction = c_i/fprime
        end
    end
end
Re_IG





