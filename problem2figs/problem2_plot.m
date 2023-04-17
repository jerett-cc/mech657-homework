close all; clear all; clc;

my_mach = csvread("problem2_mach.out");
exact_mach_b = csvread("problem2_before_shock_mach.csv")
exact_mach_a = csvread("problem2_after_shock_mach.csv")
exact_mach = [exact_mach_b; exact_mach_a];

my_pressure = csvread("problem2_pressure.out");
exact_pressure_b = csvread("problem2_before_shock_pressure.csv")
exact_pressure_a = csvread("problem2_after_shock_pressure.csv")
exact_pressure = [exact_pressure_b; exact_pressure_a];

my_density = csvread("problem2_density.out")
exact_density_b = csvread("problem2_before_shock_density.csv")
exact_density_a = csvread("problem2_after_shock_density.csv")
exact_density = [exact_density_b; exact_density_a];

my_temp = csvread("problem2_temp.out")
exact_temp_b = csvread("problem2_before_shock_temp.csv")
exact_temp_a = csvread("problem2_after_shock_temp.csv")
exact_temp = [exact_temp_b; exact_temp_a];


figure()
scatter(my_mach(:,1), my_mach(:,2), '.');
hold on
plot(exact_mach(:,1), exact_mach(:,2));
ylim([0.15, 1.6])
title("mach number, CFL = 50")
ylabel("mach")
legend(["standard implicit tmm with 90 grid points", "exact"])

figure()
scatter(my_pressure(:,1), my_pressure(:,2),'.');
hold on
plot(exact_pressure(:,1), exact_pressure(:,2), 'LineWidth' , 1);
title("pressure, CFL = 50")
ylabel("pressure in Pa")
legend(["standard implicit tmm with 90 grid points", "exact"])

figure()
scatter(my_density(:,1), my_density(:,2),'.');
hold on
plot(exact_density(:,1), exact_density(:,2), 'LineWidth' , 1);
title("density, CFL = 50")
ylabel("density (kg/m^2)")
legend(["standard implicit tmm with 90 grid points", "exact"])


figure()
scatter(my_temp(:,1), my_temp(:,2),'.');
title("temperature, CFL = 50")
hold on
plot(exact_temp(:,1), exact_temp(:,2))
ylabel("temperature (K)")
% ylim([280,300])
legend(["standard tmm with 90 gird points"])
