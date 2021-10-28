---
documentclass:
- article
geometry:
- top=1in
- left=1in
---

the following shows the matlab code for the homework assignment:

``` MATLAB
%%
 %  File: HW1.m
 % 
 %  Author: Thomas Kost
 %  
 %  Date: 21 October 2021
 %  
 %  @brief filter design hw 1
 %
 clear; clc; close all;
 
%% 3.4 All pass filter phase response
num = [0.6 1];
den = [1 0.6];
vis = fvtool(num,den);
set(vis,'Analysis','phase') 
set(vis,'Analysis','grpdelay')

%% 5.2 Impule Invariance design

w = 0.3*pi;
Omega = (2*pi)*1e7;
T = w/Omega;
num = [1-exp(-Omega*T)];
den = [1 -exp(-Omega*T)];
freqz(num,den);

w = 0.03*pi;
Omega = (2*pi)*1e7;
T = w/Omega;
num = [1-exp(-Omega*T)];
den = [1 -exp(-Omega*T)];
freqz(num,den); 
%% 6.3
 Omega_P = 48000;
 Omega_S = 80000;
 Fs = 192000;
 Omega_C = Omega_P;

 K= Omega_P/Omega_S;
 E = 0.5;
 ATT = 45;
 epsilon = sqrt((10^(E/20))^2 -1);
 A = 10^(ATT/20);
 K1 = epsilon/sqrt(A^2+1);
 N = ceil(acosh(1/K1)/acosh(1/K));
 
 alpha = 1/epsilon +sqrt(1 +1/(epsilon^2));
 a =(alpha^(1/N)-alpha^(-1/N))/2;
 b =(alpha^(1/N)+alpha^(-1/N))/2;
 
 theta = linspace(0,2*pi,2*N+1);
 theta = theta(theta>pi/2 & theta <3*pi/2);
 
 sigma = a*Omega_C*cos(theta);
 omega = b*Omega_C*sin(theta);
 poles = sigma +1j*omega;
 
 s = tf('s');
 z = tf('z', 1/Fs);
 H_c = 1;
 H_d = 1;
 T = 1/192000;
 for i=1:N
     H_d = H_d*poles(i)/(((2/T)*(z-1)/(z+1))-poles(i));
     H_c = H_c*poles(i)/(s-poles(i));
 end
 opts = bodeoptions;
 opts.freqscale = 'linear';
 bodeplot(H_c, H_d,opts)
 figure()
 freqz(cell2mat(H_d.numerator), cell2mat(H_d.denominator))
 
 
 %% 7.3 LP-HP transform
 a = 0.5*(1+sqrt(3));
 b = sqrt(3);
 num = [21 33 21];
 den = [16+4*a+b 8+17*a+8*b 1+4*a+16*b];
 tf_plt = figure();
 freqz(num,den);
 orig_plt = figure();
 freqz([1 -1 1],[1 -a b])
saveas(tf_plt, "7_3_lp_bp_transformed_plot.jpg"); 
saveas(tf_plt, "7_3_lp_bp_original_plot.jpg");

%% 8.3 Frequency response
wc1 = 0.4*pi;
wc2 = 0.5*pi;
wc = 0.25*pi;
p = cot((wc2-wc1)/2)*tan(wc/2);
lambda = cos((wc1+wc2)/2)/cos((wc1-wc2)/2);
gamma = 1-exp(-0.25*pi);
num = gamma*[p-1 -2*lambda*p p+1];
den = [p+1+exp(-0.25*pi)*(p-1) -2*lambda*p*(1+exp(-0.25*pi)) p-1+exp(-0.25*pi)*(p+1)];
lp2bp = figure();
freqz(num,den)
hold on;
freqz([1-exp(-0.25*pi)],[1 -exp(-0.25*pi)])
saveas(lp2bp, "8_3_lp_bp_plot.jpg");

 ```