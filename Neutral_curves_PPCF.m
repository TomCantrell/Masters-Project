
% =========================================================================
% -------- Plotting the neutral curve for Poisuille-Couette flow  ---------
% =========================================================================
u_W = [0.01020];%,0.0375,0.15,0.5];
% 500secs
Re_0=[];alpha_0=[];
for r=1:length(u_W)
    Re_upper=[];
    Re_lower=[];
    a_upper=[];
    a_lower=[];
    c_r_upper=[];
    c_r_lower=[];
    %alpha=1.09; % Wavenumber
    a=0.1:0.005:1.085; % Wavenumber alpha
    for p=1:length(a)
        N=59; % N+1 Chebyshev polynomials
        % Generate the approximation to the function f
        Re=3500:1000:50000;
        c_i=zeros([length(Re) 1]);
        for k=1:length(Re)
            c_i(k)= imag(spectral_PPCF(Re(k),a(p),u_W(r),N));
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
                   fprime = (imag(spectral_PPCF(Re_IG(k)+50,a(p),u_W(r),N))-imag(spectral_PPCF(Re_IG(k)-50,a(p),u_W(r),N)))/(2*50);
                   c_i = imag(spectral_PPCF(Re_IG(k),a(p),u_W(r),N));
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
            c_r_upper = [c_r_upper;real(spectral_PPCF(Re_IG(2),a(p),u_W(r),N))];
            c_r_lower = [c_r_lower;real(spectral_PPCF(Re_IG(1),a(p),u_W(r),N))];
        elseif length(Re_IG)==1
            Re_lower = [Re_lower;Re_IG(1)];
            a_lower = [a_lower;a(p)];
            c_r_lower = [c_r_lower;real(spectral_PPCF(Re_IG(1),a(p),u_W(r),N))];        
        end
        Re_IG
    end
    %% Finer part of the upper branch
    N=59;
    %Re_u = 3900:200:40000; % Reynolds numbers corresponding to the upper part of the neutral curve
    k=0;
    Re_tmp=1;
    Re_u=[];
    while Re_tmp<Re_upper(length(Re_upper)-1)
        Re_tmp = Re_lower(length(Re_lower))+k*50;
        k=k+1;
        Re_u= [Re_u;Re_tmp];
    end
    alpha_u=[];
    c_r_u = [];
    for p=1:length(Re_u)
        %alpha=0.5:0.0025:1.15;
        a_=0.2+a_lower(length(a_lower));
        alpha=[];
        while a_ > a_lower(length(a_lower)-1)
           a_ = a_-0.0025; 
           alpha=[alpha;a_]; 
        end        
        alpha=flip(alpha);
        %
        c_i_upper=zeros([length(alpha) 1]);
        for j=1:length(alpha)
            c_i_upper(j) = imag(spectral_PPCF(Re_u(p),alpha(length(alpha)-j+1),u_W(r),N));
            %alpha(length(alpha)-j+1)
        end
        for k=2:length(c_i_upper)
           if c_i_upper(k)>0 && c_i_upper(k-1)<0 
              %Re_u = [Re_u;Re];
              alpha_g = alpha(length(alpha)-k+1); % guess for alpha
              break
           end
        end
        %alpha_g;
        % =======  Newton-Raphson implementation for the upper branch  =======
        correction=10;
        while correction>10^(-4)
           % Perform iterations of Newton-Raphson root-finding algorithm until the threshold in the while loop is reached
           fprime = (imag(spectral_PPCF(Re_u(p),alpha_g+0.05,u_W(r),N))-imag(spectral_PPCF(Re_u(p),alpha_g-0.05,u_W(r),N)))/(2*0.05);
           c_i = imag(spectral_PPCF(Re_u(p),alpha_g,u_W(r),N));
           % Update estimate of root
           alpha_g = alpha_g - c_i/fprime;
           correction = c_i/fprime;
        end
        c_r_u = [c_r_u;real(spectral_PPCF(Re_u(p),alpha_g,u_W(r),N))];
        alpha_g
        Re_u(p);
        alpha_u = [alpha_u;alpha_g];
    end
    alpha_u;
    
    a_l = a_lower;
    Re_l = Re_lower;
    c_r = c_r_lower;
    endd = length(Re_u)+length(Re_upper);
    a_upper = flip(a_upper);
    Re_upper = flip(Re_upper);
    c_r_upper=flip(c_r_upper);
    for k=1:endd-1
        if k>length(Re_u)
            Re_l=[Re_l;Re_upper(k-length(Re_u)+1)];
            a_l=[a_l;a_upper(k-length(Re_u)+1)];
            c_r=[c_r;c_r_upper(k-length(Re_u)+1)];
        else
            Re_l=[Re_l;Re_u(k)];
            a_l=[a_l;alpha_u(k)];
            c_r=[c_r;c_r_u(k)]
        end
    end
    %c_r=[];
    %c_r = [c_r;c_r_lower];c_r=[c_r;c_r_u];[c_r;flip(c_r_upper)];
    % Vectors of data to be plotted
    if r==1
        Re_0 = Re_l;
        alpha_0 = a_l;
        c_r_0 = c_r;
    elseif r==2
        Re_1 = Re_l;
        alpha_1 = a_l;
        c_r_1 = c_r;
    elseif r==3
        Re_2 = Re_l;
        alpha_2 = a_l;
        c_r_2 = c_r;
    elseif r==4
        Re_3 = Re_l;
        alpha_3 = a_l;
        c_r_3 = c_r;
    elseif r==5
        Re_4 = Re_l;
        alpha_4 = a_l;
        c_r_4 = c_r;
    end
end

%% Plotting the neutral curve with the curve for wavespeed

figure(1)
plot(Re_0,alpha_0,"color","k","LineWidth", 1.45)
hold on
plot(Re_1,alpha_1,"--","color","k","LineWidth", 1.45)
hold on 
plot(Re_2,alpha_2,"-.","color","k","LineWidth", 1.45)
hold on
plot(Re_3,alpha_3,":","color","k","LineWidth", 1.45)
grid on
xlim([0 50000])
ylim([0 1.2])
xlabel("Re")
ylabel("\alpha")
lgd=legend("u_w=0.0","u_w=0.0375","u_w=0.15","u_w=0.35","Location","Best");
lgd.FontSize=10;
% figure(2)
% plot(Re_0,c_r_0,"color","k","LineWidth", 1)
% hold on
% plot(Re_1,c_r_1,"--","color","k","LineWidth", 1)
% hold on 
% plot(Re_2,c_r_2,"-.","color","k","LineWidth", 1)
% hold on
% plot(Re_3,c_r_3,":","color","k","LineWidth", 1)
% grid on
% xlim([0 40000])
% ylim([0 0.5])
% xlabel("Re")
% ylabel("c_r")
% legend("u_w=0.0","u_w=0.0375","u_w=0.15","u_w=0.25","Location","Best")

















