close all; clear all; clc;

C = csvread("problem3_mach.out");
C1 = csvread("velocity_number.txt");



my_pressure = csvread("problem3_pressure.out");
exact_pressure = csvread("pressure_number.txt")


figure()
plot(Bnew(:,1), C(:,1));
hold on
plot(C1(:,2), C1(:,1));
ylabel("velocity")
legend(["20 grid solution", "actual"])

figure()
% plot(my_pressure(:,1), my_pressure(:,2), 'linewidth',2);
scatter(my_pressure(:,1), my_pressure(:,2),'.');
hold on
% plot(exact_pressure(:,2), exact_pressure(:,1), 'LineWidth' , 1);
ylabel("pressure in Pa")
legend(["Computed solution 99 grid points", "Exact solution"])



% C = csvread("step-99_temp.out");
% 
% figure()
% plot(Bnew(:,1), C(:,1));
% ylabel("temperature in K")

D = csvread("problem3_density.out");
D1 = csvread("density_number.txt");

figure()
plot(D(:,1), D(:,2));
hold on
% plot(D1(:,2), D1(:,1));
% legend(["calculated", "actual"])
ylabel("density (kg/m^3)")