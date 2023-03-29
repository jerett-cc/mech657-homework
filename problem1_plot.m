close all; clear all; clc;

C = csvread("step-4_velocity.out");
C1 = csvread("velocity_number.txt");



B = csvread("step-4_pressure.out");
Bnew = [ 0, 97239.4
 0.47619, 96514.9
 0.952381, 95585.7
 1.42857, 94404.1
 1.90476, 92925.2
 2.38095, 91122
 2.85714, 89012.4
 3.33333, 86696.8
 3.80952, 84397.6
 4.28571, 82473.4
 4.7619, 81353.5
 5.2381, 81254.6
 5.71429, 81645
 6.19048, 82384.6
 6.66667, 83402.3
 7.14286, 84613.3
 7.61905, 85934
 8.09524, 87291.8
 8.57143, 88629.2
 9.04762, 89904.9
 9.52381, 91092.1
 10, 92176]

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

figure()
plot(C(:,1), C(:,2));
hold on
plot(C1(:,2), C1(:,1));
ylabel("velocity")
legend(["20 grid solution", "actual"])

figure()
plot(B(:,1), B(:,2));
hold on
plot(Bnew(2:21,1), P);
ylabel("pressure in Pa")
legend(["actual", "20 grid solution"])



% C = csvread("step-99_temp.out");
% 
% figure()
% plot(Bnew(:,1), C(:,1));
% ylabel("temperature in K")

D = csvread("step-4_density.out");
D1 = csvread("density_number.txt");

figure()
plot(D(:,1), D(:,2));
hold on
plot(D1(:,2), D1(:,1));
legend(["calculated", "actual"])
ylabel("density (kg/m^3)")