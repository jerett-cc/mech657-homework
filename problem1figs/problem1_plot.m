close all; clear all; clc;

my_mach = csvread("problem1_mach.out");
exact_mach = csvread("problem1_mach.csv");

my_pressure = csvread("problem1_pressure.out");
exact_pressure = csvread("problem1_pressure.csv")

my_density = csvread("problem1_density.out")
exact_density = csvread("problem1_density.csv")

my_temp = csvread("problem1_temp.out")
exact_temp = csvread("problem1_temp.csv")

my_vel = csvread("problem1_velocity.out")
exact_vel = csvread("problem1_velocity.csv")

%initial
my_machINIT = csvread("problem1INIT_mach.out");
exact_mach = csvread("problem1_mach.csv");

my_pressureINIT = csvread("problem1INIT_pressure.out");
exact_pressure = csvread("problem1_pressure.csv")

my_densityINIT = csvread("problem1INIT_density.out")
exact_density = csvread("problem1_density.csv")

my_tempINIT = csvread("problem1INIT_temp.out")
exact_temp = csvread("problem1_temp.csv")

my_velINIT = csvread("problem1INIT_velocity.out")
exact_vel = csvread("problem1_velocity.csv")
%end initial


figure()
scatter(my_mach(:,1), my_mach(:,2), '.');
hold on
plot(exact_mach(:,1), exact_mach(:,2));
ylim([0.15,0.7])
title("mach number, CFL = 100")
ylabel("mach")
legend(["standard implicit tmm with 99 grid points", "exact"])

figure()
scatter(my_pressure(:,1), my_pressure(:,2),'.');
hold on
plot(exact_pressure(:,1), exact_pressure(:,2), 'LineWidth' , 1);
title("pressure, CFL = 100")
ylabel("pressure in Pa")
legend(["standard implicit tmm with 99 grid points", "exact"])

figure()
scatter(my_density(:,1), my_density(:,2),'.');
hold on
plot(exact_density(:,1), exact_density(:,2), 'LineWidth' , 1);
title("density, CFL = 100")
ylabel("density (kg/m^2)")
legend(["standard implicit tmm with 99 grid points", "exact"])


figure()
scatter(my_temp(:,1), my_temp(:,2),'.');
title("temperature, CFL = 100")
hold on
plot(exact_temp(:,1), exact_temp(:,2))
ylabel("temperature (K)")
ylim([280,300])
legend(["standard tmm with 99 grid points"])

figure()
scatter(my_vel(:,1), my_vel(:,2),'.');
title("velocity")
hold on
plot(exact_vel(:,1), exact_vel(:,2))
ylabel("velocity (m/s)")
% ylim([280,300])
legend(["standard tmm with 99 grid points"])

test_err2 = zeros(1,7);
test_iter2 = zeros(1,7);
test_err2(1) = 65000;
test_iter2(1) = 1;
for iii=2:7
    test_err2(iii) = test_err2(iii-1)/4;
    test_iter2(iii) = test_iter2(iii-1)*2;
end

% figure()
% %semilogy(error)
% hold on
% semilogy(error100)
% semilogy(error1000)
% % loglog(test_iter2 , test_err2);
% title("error")
% ylabel("l2 error")
% xlabel("iterations")
% legend(["CFL=50", "CFL=100", "CFL=1000"])

%initial
figure()
scatter(my_machINIT(:,1), my_machINIT(:,2), '.');
hold on
plot(exact_mach(:,1), exact_mach(:,2));
ylim([0.15,0.7])
title("mach")

figure()
scatter(my_pressureINIT(:,1), my_pressureINIT(:,2),'.');
hold on
plot(exact_pressure(:,1), exact_pressure(:,2), 'LineWidth' , 1);
title("press")

figure()
scatter(my_densityINIT(:,1), my_densityINIT(:,2),'.');
hold on
plot(exact_density(:,1), exact_density(:,2), 'LineWidth' , 1);
title("density")


figure()
scatter(my_tempINIT(:,1), my_tempINIT(:,2),'.');
title("temperature")
hold on
plot(exact_temp(:,1), exact_temp(:,2))

figure()
scatter(my_velINIT(:,1), my_velINIT(:,2),'.');
title("velocity")
hold on
plot(exact_vel(:,1), exact_vel(:,2))
