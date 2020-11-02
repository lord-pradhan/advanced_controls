clear 
clc
close all

load('lab3_ws.mat')

%% init variables
m_c = 0.493; m_p = 0.312; I_p = 0.00024;
l = 0.04; f = 0.01; k_t = 0.11; R = 10; r = 0.0335; g=9.81;
sample_time = 0.005;

%% state space
X_ref = [0,0,0.0,0]';

A = [0,1,0,0;
    0, -(I_p + m_p*l^2)*f/( I_p*(m_c+m_p) + m_c*m_p*l^2 ), m_p^2 * g* l^2/(I_p*(m_c+m_p) + m_c*m_p*l^2), 0;
    0,0,0,1;
    0, -m_p*l*f/( I_p*(m_c+m_p) + m_c*m_p*l^2 ), m_p*g*l*(m_c + m_p)/( I_p*(m_c+m_p) + m_c*m_p*l^2 ), 0];
B = [0;
    (I_p + m_p*l^2)/( I_p*(m_c+m_p) + m_c*m_p*l^2 );
    0;
    m_p*l/( I_p*(m_c+m_p) + m_c*m_p*l^2 )] * 2*k_t/(R*r);

% C = [1,0,0,0;
%     0,0,1,0]; 
C = eye(4);

D = zeros(size(C,1), size(B,2));

[n,d] = ss2tf(A,B,C,D);

%% lqr
Q = diag([0.001, 0.1, 5000, 1]);
R = 0.2;

[k,S,CLP] = lqr(A,B,Q,R);
k(4) = 4.0; % reduce gains slightly to account for gyro noise

%gains that work - {1.03, -113, -6.2}, {1.03, -154.6, -3}

%% closed loop TF
[n_closed, d_closed] = ss2tf(A - B*k, B, C, D);

cltf_X = tf(n_closed(1,:), d_closed);
cltf_theta = tf(n_closed(3,:), d_closed);


%% 
% time_real = 0:0.005:2.5;
% figure(1)
% xlabel('Time (s)')
% yyaxis left
% plot(time_real, freq4twice(1:501,1))
% set(gca, 'YLim', [min(freq4twice(1:501,1))*1.1, max(freq4twice(1:501,1))*1.1]);
% ylabel(gca, 'Forcing input (V)')
% 
% yyaxis right
% plot(time_real, freq4twice(1:501, 3), 'r-')
% set(gca, 'YLim', [min(freq4twice(1:501,3)), max(freq4twice(1:501,3))]);
% ylabel(gca, 'Position output (m)')
% 
% 
% figure(2)
% xlabel('Time (s)')
% yyaxis left
% plot(time_real, freq4twice(1:501,1))
% set(gca, 'YLim', [min(freq4twice(1:501,1))*1.1, max(freq4twice(1:501,1))*1.1]);
% ylabel(gca, 'Forcing input (V)')
% 
% yyaxis right
% plot(time_real, freq4twice(1:501, 4), 'r-')
% set(gca, 'YLim', [min(freq4twice(1:501,4))*1.1, max(freq4twice(1:501,4))*1.1]);
% ylabel(gca, 'Angle output (rad)')

%% bode plot
amp_gain_pos = [0.035/8, 0.03/8];
amp_gain_angle = [0.042/8, 0.05/8];
% 
phase_angle = [pi/2+0.35, pi/2+0.5] - pi/2*ones(5,1);
phase_pos = [0, 5]*pi/180;


freqs = 2*pi*[1.687, 3];
Fs = 66.67;
L =  2401;
gainX_data = zeros(1,length(freqs));
gainTheta_data = zeros(1,length(freqs));
phaseTheta_data = zeros(1,length(freqs));
phaseX_data = zeros(1,length(freqs));


logs = zeros(L, 4, length(freqs));
logs(:,:,1) = freq3twice(1:L,:); 
logs(:,:,2) = freq4twice(1:L,:);

logs(isnan(logs))=0;

for i = 1:length(freqs)
    log_temp = logs(:,:,i);
    fid = round( freqs(i)*L/(2*pi*Fs))+1;
    U = goertzel(-log_temp(1:L,1)', fid);
    pos = goertzel(log_temp(1:L,3)', fid);
    angleTemp = goertzel(log_temp(1:L, 4)', fid);
    
    gainX_data(i) = abs(pos/U);
    gainTheta_data(i) = abs(angleTemp/U);
    phaseX_data(i) = angle(pos/U);
    phaseTheta_data(i) = angle(angleTemp/U);
end

fin_gainX = gainX_data; %amp_gain_pos
fin_gainTheta = gainTheta_data; %amp_gain_angle
fin_phaseX = phaseX_data; %phase_pos
fin_phaseTheta = phaseTheta_data; %phase_angle
beginPos = 1;

%% plots
[magX,phaseX] = bode(cltf_X, {freqs(1)-1, freqs(end)+1});
magX = squeeze(magX); phaseX = squeeze(phaseX);
wradX = linspace(freqs(1)-1, freqs(end)+1, length(magX));

figure(1)
subplot(2,1,1)
semilogx(wradX, magX)
ylabel('Amplitude gain (ratio)')
title('Sim vs real bode plot for position')
hold on
plot(freqs(beginPos:end),fin_gainX, 'r*', freqs, amp_gain_pos, 'kx')
legend('TF', 'Goertzel', 'Inspection')
hold off
grid
subplot(2,1,2)
semilogx(wradX, phaseX)
xlabel('Frequency (rad/s)')
ylabel('Phase lag (deg)')
hold on
plot(freqs(beginPos:end), 180/pi*fin_phaseX, 'r*', freqs, phase_pos, 'kx')
hold off
grid

[magTheta,phaseTheta] = bode(cltf_theta, {freqs(1)-1, freqs(end)+1});
magTheta = squeeze(magTheta); phaseTheta = squeeze(phaseTheta);
wradTheta = linspace(freqs(1)-1, freqs(end)+1, length(magTheta));

figure(2)
subplot(2,1,1)
semilogx(wradTheta, magTheta)
ylabel('Amplitude gain (ratio)')
title('Sim vs real bode plot for angle')
hold on
plot(freqs,fin_gainTheta, 'r*', freqs, amp_gain_angle, 'kx')
legend('TF', 'Goertzel', 'Inspection')
hold off
grid
subplot(2,1,2)
semilogx(wradTheta, phaseTheta)
xlabel('Frequency (rad/s)')
ylabel('Phase lag (deg)')
hold on
plot(freqs, 180/pi*fin_phaseTheta, 'r*', freqs, phase_angle, 'kx')
hold off
grid

%%
%% Simulation
% amplitude = 01;
% 
% t = 0:0.01:5;
% u = amplitude*sin(freqs(5)*t);
% Y_pos = lsim(cltf_X, u, t);
% Y_theta = lsim(cltf_theta, u, t);

% figure(3)
% plot(t, Y_pos)
% title('Closed loop position sim')
% 
% figure(4)
% plot(t, Y_theta)
% title('Closed loop angle sim')

