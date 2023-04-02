close all; clear all; clc;

C = csvread("problem1_velocity.out");
C1 = csvread("velocity_number.txt");



my_pressure = csvread("problem1_pressure.out");
exact_pressure = csvread("pressure_number.txt")

P =[
99319.88599805532
94916.18950579465
96581.67986116301
92309.78341962787
91766.89133877688
88498.53867728838
86140.06085491576
83338.91175996303
81011.82840235543
79678.31102781248
79878.87469736804
80566.27882523059
81030.03624773549
82635.98650444369
83416.07973278993
85631.24552573537
86395.58966309673
88937.01792202386
89704.74740239118
92227.95877175649
    ]

% figure()
% plot(Bnew(:,1), C(:,1));
% hold on
% plot(C1(:,2), C1(:,1));
% ylabel("velocity")
% legend(["20 grid solution", "actual"])

figure()
plot(my_pressure(:,1), my_pressure(:,2));
hold on
plot(exact_pressure(:,2), exact_pressure(:,1));
ylabel("pressure in Pa")
legend(["mine", "20 grid solution"])



% C = csvread("step-99_temp.out");
% 
% figure()
% plot(Bnew(:,1), C(:,1));
% ylabel("temperature in K")

D = csvread("problem1_density.out");
D1 = csvread("density_number.txt");

figure()
plot(D(:,1), D(:,2));
hold on
plot(D1(:,2), D1(:,1));
legend(["calculated", "actual"])
ylabel("density (kg/m^3)")