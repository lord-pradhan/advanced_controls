%Author - Roshan Pradhan
% This script is written for lab3 of 24774 course at CMU

clear 
clc
close all

load('lab3_ws.mat')
%% Init variables
m_c = 0.493; m_p = 0.312; I_p = 0.00024;
l = 0.04; f = 0.01; k_t = 0.11; R = 10; r = 0.0335; g=9.81;
sample_time = 0.005;

%% Nominal state space model
X_ref = [0,0,0.0,0]';

A = [0,1,0,0;
    0, -(I_p + m_p*l^2)*f/( I_p*(m_c+m_p) + m_c*m_p*l^2 ),...
    m_p^2 * g* l^2/(I_p*(m_c+m_p) + m_c*m_p*l^2), 0;
    0,0,0,1;
    0, -m_p*l*f/( I_p*(m_c+m_p) + m_c*m_p*l^2 ), ...
    m_p*g*l*(m_c + m_p)/( I_p*(m_c+m_p) + m_c*m_p*l^2 ), 0];
B = [0;
    (I_p + m_p*l^2)/( I_p*(m_c+m_p) + m_c*m_p*l^2 );
    0;
    m_p*l/( I_p*(m_c+m_p) + m_c*m_p*l^2 )] * 2*k_t/(R*r);
C = eye(4);
D = zeros(size(C,1), size(B,2));

[n,d] = ss2tf(A,B,C,D);

%% LQR design
Q = diag([0.001, 0.1, 5000, 1]);
R = 0.2;

[k,S,CLP] = lqr(A,B,Q,R);
k(4) = 4.0; % reduce gains slightly to account for gyro noise

%% Closed loop TF
[n_closed, d_closed] = ss2tf(A - B*k, B, C, D);

cltf_X = tf(n_closed(1,:), d_closed);
cltf_theta = tf(n_closed(3,:), d_closed);

figure(1)
bode(cltf_X)
title('Bode plot for position CLTF')

figure(2)
bode(cltf_theta)
title('Bode plot for angle CLTF')

%% Process freq data
time_real = 0:0.005:2;
figure(3)
xlabel('Time (s)')
yyaxis left
plot(time_real, -freq4twice(1:401,1))
set(gca, 'YLim', [min(freq4twice(1:401,1))*1.1, max(freq4twice(1:401,1))*1.1]);
ylabel(gca, 'Forcing input (V)')

yyaxis right
plot(time_real, freq4twice(1:401, 3), 'r-')
set(gca, 'YLim', [min(freq4twice(1:401,3)), max(freq4twice(1:401,3))]);
ylabel(gca, 'Position output (m)')


figure(4)
xlabel('Time (s)')
yyaxis left
plot(time_real, -freq4twice(1:401,1))
set(gca, 'YLim', [min(freq4twice(1:401,1))*1.1, max(freq4twice(1:401,1))*1.1]);
ylabel(gca, 'Forcing input (V)')

yyaxis right
plot(time_real, freq4twice(1:401, 4), 'r-')
set(gca, 'YLim', [min(freq4twice(1:401,4))*1.1, max(freq4twice(1:401,4))*1.1]);
ylabel(gca, 'Angle output (rad)')

%% Frequency domain - Visual inspection vs Goertzel
freqs = 2*pi* logspace(log10(0.3), log10(3), 5); % in rad/s

% visual inspected values
amp_gain_angle = [0.027/4, 0.025/4, 0.032/4, 0.025/4, 0.022/4];
amp_gain_pos = [0.035/4, 0.028/4, 0.04/4, 0.03/2, 0.01/2];
phase_angle = [pi/2, pi/2, pi/2+0.1745, pi/2+0.35, pi/2+0.5] - pi/2*ones(5,1);
phase_pos = [0, 0, -5, 0, 5]*pi/180;

% Goerztel script
Fs = 66.67;
L =  1401;
logs = zeros(L, 4, length(freqs));
logs(:,:,1) = freq0(1:L,:); 
logs(:,:,2) = freq1(1:L,:);
logs(:,:,3) = freq2(1:L,:);
logs(:,:,4) = freq3(1:L,:); 
logs(:,:,5) = freq4(1:L,:);

fin_gainX = zeros(1,length(freqs));
fin_gainTheta  = zeros(1,length(freqs));
fin_phaseTheta = zeros(1,length(freqs));
fin_phaseX = zeros(1,length(freqs));

for i = 1:length(freqs)
    log_temp = logs(:,:,i);
    fid = round( freqs(i)*L/(2*pi*Fs))+1;
    U = goertzel(-log_temp(1:L,1)', fid);
    pos = goertzel(log_temp(1:L,3)', fid);
    angleTemp = goertzel(log_temp(1:L,4)', fid);
    
    fin_gainX(i) = abs(pos/U);
    fin_gainTheta(i) = abs(angleTemp/U);
    fin_phaseX(i) = angle(pos/U);
    fin_phaseTheta(i) = angle(angleTemp/U);
end


%% Plots for Bode comparison

% nominal model
[magX,phaseX] = bode(cltf_X, {freqs(1)-1, freqs(end)+1});
magX = squeeze(magX); phaseX = squeeze(phaseX);
wradX = linspace(freqs(1)-1, freqs(end)+1, length(magX));

%identified model
[sys_magX,sys_phaseX] = bode(tfPosFin, {freqs(1)-1, freqs(end)+1});
sys_magX = squeeze(sys_magX); sys_phaseX = squeeze(sys_phaseX);
sys_wradX = linspace(freqs(1)-1, freqs(end)+1, length(sys_magX));


figure(5)
subplot(2,1,1)
semilogx(wradX, magX)
ylabel('Amplitude gain (ratio)')
title('Sim vs real bode plot for position')
hold on
plot(freqs(1:end),fin_gainX, 'r*', freqs, amp_gain_pos, 'kx')
semilogx(sys_wradX, sys_magX, 'k--')
legend('TF', 'Goertzel', 'Inspection', 'Identified model')
hold off
grid
subplot(2,1,2)
semilogx(wradX, phaseX)
xlabel('Frequency (rad/s)')
ylabel('Phase lag (deg)')
hold on
plot(freqs(1:end), 180/pi*fin_phaseX, 'r*', freqs, phase_pos, 'kx')
semilogx(sys_wradX, sys_phaseX, 'k--')
hold off
grid

% nominal model
[magTheta,phaseTheta] = bode(cltf_theta, {freqs(1)-1, freqs(end)+1});
magTheta = squeeze(magTheta); phaseTheta = squeeze(phaseTheta);
wradTheta = linspace(freqs(1)-1, freqs(end)+1, length(magTheta));

%identified model
[sys_magTheta,sys_phaseTheta] = bode(tfAngleFin, {freqs(1)-1, freqs(end)+1});
sys_magTheta = squeeze(sys_magTheta); sys_phaseTheta = squeeze(sys_phaseTheta);
sys_wradTheta = linspace(freqs(1)-1, freqs(end)+1, length(sys_magTheta));


figure(6)
subplot(2,1,1)
semilogx(wradTheta, magTheta)
ylabel('Amplitude gain (ratio)')
title('Sim vs real bode plot for angle')
hold on
plot(freqs,fin_gainTheta, 'r*',  freqs, amp_gain_angle, 'kx')
semilogx(sys_wradTheta, sys_magTheta, 'k--')
legend('TF', 'Goertzel', 'Inspection', 'Identified model')
hold off
grid
subplot(2,1,2)
semilogx(wradTheta, phaseTheta)
xlabel('Frequency (rad/s)')
ylabel('Phase lag (deg)')
hold on
plot(freqs, 180/pi*fin_phaseTheta, 'r*',  freqs, phase_angle, 'kx')
semilogx(sys_wradTheta, sys_phaseTheta, 'k--')
hold off
grid


