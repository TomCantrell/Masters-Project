
% Plot of the first 5 Chebyshev polynomials

z = -1:0.025:1; 
T_0 = chebyshevT(0,z);
T_1 = chebyshevT(1,z);
T_2 = chebyshevT(2,z);
T_3 = chebyshevT(3,z);
T_4 = chebyshevT(4,z);
T_5 = chebyshevT(5,z);
% Plotting
%plot(z,T_4,"-o",z,T_3,"-.",z,T_2,"--",z,T_1,"*-",z,T_0)
plot(z,T_5,"-x","LineWidth",1.0,"color",[0 0.4470 0.7410])
hold on
plot(z,T_4,"-o","LineWidth",1.0,"color",[0.3010 0.7450 0.9330])
hold on
plot(z,T_3,"-.","LineWidth",1.0,"color",[0.4660 0.6740 0.1880])

hold on
plot(z,T_2,"*-","LineWidth",1.0,"color",[0 0 1])

hold on
plot(z,T_1,"--","LineWidth",1.0,"color","m")

hold on

plot(z,T_0,"LineWidth",1.0,"color",[0 0 0])
legend("T_5(z)","T_4(z)","T_3(z)","T_2(z)","T_1(z)","T_0(z)","Location","Best")
xlim([-1.2 1.2])
ylim([-1.2 1.2])
grid on
xlabel("z")
grid on



