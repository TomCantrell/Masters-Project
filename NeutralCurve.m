
% =======================================================
% --------  Neutral curve of stability - PPF    ---------
% =======================================================

% Computation of the neutral curve of linear stability obtained by
% utilising the local iterative method derived from the finite-difference
% implementation. 

% ===================== Inital Guess ====================
c_r = [];
Re = 5772.22;
alpha = 1.02056;
m=4;
n = 2^m;
t = 2^(m-1); % # grid points
h = 2/(t-1);  
U = zeros([t 1]);
for i = 1:t-1
    U(i)= 1-(-1 + h*(i-1))^2;
end
A = zeros(n);
B = zeros(n);
b = sqrt(-1)/(alpha*Re);
% phi' = 0 BC
A(1,2) = -3/(2*h);
A(1,4) = 4/(2*h);
A(1,6) = -1/(2*h);
% phi = 0 BC
A(2,2) = 1;
for j = 2:n/2 - 1
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
%------------------------------------------------------------------------
% interpolation of eigenvector to approximate to a larger grid
% interpolate until size of Psiginter exceeds 
o=1;
Psiginter = 1;
while length(Psiginter) < 1000
    Psiginter = zeros([2*(length(Psig)/2 -1)+n 1]);
    i=1;
    for k=1:4:length(Psiginter)    
        Psiginter(k)=  Psig(i); % Psi eqn
        Psiginter(k+1)=  Psig(i+1); % Phi eqn
        i=i+2;
    end
    p=1;
    for j=3:4:length(Psiginter)
        Psiginter(j)= (Psig(p)+Psig(p+2))/2; %
        Psiginter(j+1)= (Psig(p+1)+Psig(p+3))/2; %
        p=p+2;
    end
    Psig = Psiginter;
    n = length(Psiginter);
    o=o+1;
end
N = length(Psiginter);
PsiG = zeros([N+1 1]);
PsiG(1) = cG;
for j=2:N+1
    PsiG(j) = Psiginter(j-1);
end
% To create the plot demonstrating the relationship between wave speed (the
% real part of the eigenvalue) and Re# we need to note down the wave speeds
% at different Re#s
c_r_upper = [];
c_r_lower = [];  % Initialising vectors which will hold wave speed values


% =========================================================================
InitialPsi = PsiG;
Relim = 40000; 
% ---------------------------  lower branch  ------------------------------
% Finer mesh
% Start at Re=6000, alpha=0.9, close to the lower branch but outside it in
% the stable region, step up to the curve of neutral stability. Then begin
% tracking "backwards" in terms of Re# as well as "upwards" alpha.
Re = 6000;
alpha = 0.9;
c = PsiG(1);
while abs(imag(c))>0.0005
    e = Newton(Re,alpha,PsiG);
    c = e(1);
    alpha = alpha+0.001;
    PsiG = e;
end 
PsiG_lower = PsiG;
% Tracking backwards means have to initialise a new vector and then reverse
% the order when adding it to main Re_lower vector
Re_lower_fine = [];
alpha_lower_fine = [];
c_r_lower_fine = [];
var = -1;
while var < 0
    e = Newton(Re,alpha,PsiG);
    Re;
    c_i = imag(e(1));
    if abs(c_i) < 0.0005
        Re_lower_fine = [Re_lower_fine; Re];
        alpha_lower_fine = [alpha_lower_fine; alpha];
        Re = Re-20;
        c_r_lower_fine = [c_r_lower_fine;real(e(1))];
    else
        w = Newton(Re,alpha-0.0005,PsiG);
        alpha = alpha + 0.0005
        if c_i - imag(w(1)) < 0
            var = 1;
        end
    end
    PsiG = e;
end
RElower=[];
Alphalower=[];
c_r_lower=[]
% Add the finer part of the lower branch to the full lower branch vector,
% ie just reverse the order and add to the vector
for i=1:length(Re_lower_fine)
    len = length(Re_lower_fine);
    RElower(i) = Re_lower_fine(len-i+1);
    Alphalower(i) = alpha_lower_fine(len-i+1);
    c_r_lower(i) = c_r_lower_fine(len-i+1);
end

% Initialising for remainder of lower branch
Re = RElower(length(Re_lower_fine))
alpha = Alphalower(length(alpha_lower_fine))
PsiG = Newton(Re,alpha,PsiG);
while Re < Relim  % Original code, unaltered
    v = Newton(Re,alpha,PsiG);
    Re
    if abs(imag(v(1)))<0.0005
        %Add Re and alpha to a list, thereby starting a list of alphas and
        %Re numbers where c_i =~ 0.
        RElower(length(RElower)+1)= Re; %[RElower; Re]; % NEED TO CHANGE THESE CONCATENATION DOESNT WORK ls(end+1) = element;
        Alphalower(length(Alphalower)+1) = alpha;   %[Alphalower; alpha];
        alpha = alpha - 0.001;
        % add to c_r vector
        c_r_lower(length(c_r_lower)+1) = real(v(1)); % [c_r_lower; real(v(1))];
    else
        Re = Re + 25;
    end
    PsiG = v;
end

% --------------------------  Upper branch  -------------------------------
PsiG= InitialPsi;
Re=RElower(1); % Start from left most (alpha,Re) combination
alpha=Alphalower(1);
REupper=[];
Alphaupper=[];
x=-1;
while x<0
    v=Newton(Re,alpha,PsiG);
    c_i = imag(v(1));
    if abs(c_i) < 0.0005
        REupper=[REupper;Re];
        Alphaupper=[Alphaupper;alpha];
        alpha=alpha+0.0005;
        % add to c_r vector
        c_r_upper = [c_r_upper;real(v(1))];
    else
        w = Newton(Re-20,alpha,PsiG);
        Re=Re+20;
        if c_i - imag(w(1)) < 0
            x=1;
        end
    end
    PsiG = v;
end

Re = REupper(length(REupper));
alpha = Alphaupper(length(Alphaupper));
PsiG = v;
% FINER MESH ON PART OF UPPER BRANCH
while Re<Relim
    w = Newton(Re,alpha,PsiG);
    if imag(w(1)) > 0
        Re=Re+100;%+250
    end
    if abs(imag(w(1)))<0.0005 && imag(w(1))<0
        REupper=[REupper;Re];
        Alphaupper=[Alphaupper;alpha];
        Re=Re+100; %+250
        c_r_upper = [c_r_upper;real(w(1))];
%         alpha=alpha-0.001;
    else
        alpha=alpha-0.0005;
    end
    PSiG=w;
    PsiG=PSiG;
end
%% Plotting the curves 
% ====================  Plotting the Neutral Curve  =======================
figure(1)
plot(REupper,Alphaupper,"-",RElower,Alphalower,"-", "color","k","LineWidth",1.1)
%title("Neutral Stability Curve")
xlabel("Re")
xlim([0 Relim])
ylim([0 1.2])
ylabel("\alpha")
grid on 
%text(10000,0.95, "Region of instability", "color","k")
%text(4000,0.6, "Stable", "color","k")

%% =====================   Dispersion curves   =============================
figure(2)



for Reynold = 5000:2000:11000
   c_i = [];
   PsiG = Newton(Reynold,0.25,PsiG);
   for alpha=0.25:0.05:1.35
       PsiG = Newton(Reynold,alpha,PsiG);
       c_i = [c_i;imag(PsiG(1))];    
   end
   alpha=0.25:0.05:1.35;
   if Reynold==11000 % Not equal
        plot(alpha, c_i,"-","Linewidth",0.8)
        hold on
   else
       plot(alpha, c_i,"-.","LineWidth",0.8)
       hold on
   end
end
xlim([0.2 1.4])
xL = xlim;
yL = ylim;
hold on
line(xL, [0 0],"color","k"); % x-axis
line([0.2 0.2], yL,"color","k"); % y-axis
grid on
xlabel("\alpha")
ylabel("c_i")
legend("Re=5000","Re=7000","Re=9000","Re=11000", "Location","Best")


%% ========================================================================
% Wave speed before the onset of instability 

% =========  Plotting the wave speed with the Reynolds number  ============
figure(3)
plot(REupper,c_r_upper,"-",RElower,c_r_lower,"-","color","k","LineWidth",0.9)
xlabel("Re")
ylabel("c_r")
xlim([0 Relim])
ylim([0 0.3])
grid on

% ~ 6 mins
maximum = 0;
for a = 1:length(Alphaupper)
   if  Alphaupper(a) > maximum
       a
       maximum = Alphaupper(a);
   end
end
max = 0;
for d=1:length(c_r_upper)
   if c_r_upper(d)>max
       d
       max = c_r_upper(d);
   end
end
% 102nd entry of c_r_upper


% "New" method
% ===================== Inital Guess ====================
% c_r = [];
% Re = 5772.22;
% alpha = 1.02056;
% m=7;
% n = 2^m;
% t = 2^(m-1); % # grid points
% h = 2/(t-1);  
% U = zeros([t 1]);
% for i = 1:t-1
%     U(i)= 1-(-1 + h*(i-1))^2;
% end
% A = zeros(n);
% B = zeros(n);
% b = sqrt(-1)/(alpha*Re);
% % phi' = 0 BC
% A(1,2) = -3/(2*h);
% A(1,4) = 4/(2*h);
% A(1,6) = -1/(2*h);
% % phi = 0 BC
% A(2,2) = 1;
% for j = 2:n/2 - 1
%     k = 2*j-1;
%     B(k,k) = 1;
%     A(k,k-2) = b/h^2 ; %psi_j-1
%     A(k,k) = U(j) - (2/h^2 + alpha^2)*b ; %psi_j
%     A(k,k+1) = 2 ; %phi_j
%     A(k,k+2) = b/h^2; %psi_j+1
%     k = k+1;
%     A(k,k-2) = -1/h^2;
%     A(k,k-1) = 1;
%     A(k,k) = 2/h^2 + alpha^2; 
%     A(k,k+2) = -1/h^2;
% end
% % phi' = 0 BC
% A(n-1,n-4) = 1/(2*h);
% A(n-1,n-2) = -4/(2*h);
% A(n-1,n) = 3/(2*h);
% % phi = 0 BC
% A(n,n) = 1;
% [V,D] = eig(A,B);
% eigsort = sort(diag(D));
% v = eigsort;
% cG = v(1);  
% for i=1:n
%     if D(i,i) == cG
%         k=i;
%     end
% end
% Psig = V(:,k:k); % Obtains eigenvector for corresponding eigenvalue
% %------------------------------------------------------------------------
% % interpolation of eigenvector to approximate to a finer grid
% % interpolate until size of Psiginter exceeds 
% o=1;
% Psiginter = 1;
% Psig_interpolated = 1;
% while length(Psig_interpolated) < 2^12
%     Psiginter = zeros([2*(length(Psig)/2 -1)+n 1]);
%     i=1;
%     for k=1:4:length(Psiginter)    
%         Psiginter(k)=  Psig(i); % Psi eqn
%         Psiginter(k+1)=  Psig(i+1); % Phi eqn
%         i=i+2;
%     end
%     p=1;
%     for j=3:4:length(Psiginter)
%         Psiginter(j)= (Psig(p)+Psig(p+2))/2; %
%         Psiginter(j+1)= (Psig(p+1)+Psig(p+3))/2; %
%         p=p+2;
%     end
%     
%     r = log(length(Psig))/log(2);
%     Psig_interpolated = zeros([2^(r+1) 1]);
%     for var1=1:length(Psig_interpolated)/2-1
%         Psig_interpolated(var1) = Psiginter(var1);
%     end
%     Psig_interpolated(length(Psig_interpolated)/2) = Psiginter(var1+1);
%     Psig_interpolated(length(Psig_interpolated)/2+1) = Psiginter(var1+2);
%     for var2 = (length(Psig_interpolated)/2)+2:length(Psig_interpolated)
%         Psig_interpolated(var2) = Psiginter(var2-2);
%     end
%     
%     Psig = Psig_interpolated;
%     n = length(Psig);
%     %o=o+1;
% end
% 
% N = length(Psig_interpolated);
% PsiG = zeros([N+1 1]);
% PsiG(1) = cG;
% for j=2:N+1
%     PsiG(j) = Psig_interpolated(j-1);
% end






