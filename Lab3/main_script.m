clear 
clc
close all

load('freqs_data.mat')
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

%gains that work - {1.03, -113, -6.2}, {1.03, -154.6, -3}

%% closed loop TF
[n_closed, d_closed] = ss2tf(A - B*k, B, C, D);

cltf_X = tf(n_closed(1,:), d_closed);
cltf_theta = tf(n_closed(3,:), d_closed);

%%
freqs = 2*pi* logspace(log10(0.3), log10(3), 5); % in rad/s
amplitude = 01;

t = 0:0.01:5;
u = amplitude*sin(freqs(5)*t);
Y_pos = lsim(cltf_X, u, t);
Y_theta = lsim(cltf_theta, u, t);

% figure(3)
% plot(t, Y_pos)
% title('Closed loop position sim')
% 
% figure(4)
% plot(t, Y_theta)
% title('Closed loop angle sim')

%% Process freq data
% time_real = 0:0.005:8;
% figure(1)
% xlabel('Time (s)')
% yyaxis left
% plot(time_real, freq0(1:1601,1))
% set(gca, 'YLim', [min(freq0(1:1601,1))*1.1, max(freq0(1:1601,1))*1.1]);
% ylabel(gca, 'Forcing input (V)')
% 
% yyaxis right
% plot(time_real, freq0(1:1601, 3), 'r-')
% set(gca, 'YLim', [min(freq0(1:1601,3)), max(freq0(1:1601,3))]);
% ylabel(gca, 'Position output (m)')
% 
% 
% figure(2)
% xlabel('Time (s)')
% yyaxis left
% plot(time_real, freq0(1:1601,1))
% set(gca, 'YLim', [min(freq0(1:1601,1))*1.1, max(freq0(1:1601,1))*1.1]);
% ylabel(gca, 'Forcing input (V)')
% 
% yyaxis right
% plot(time_real, freq0(1:1601, 4), 'r-')
% set(gca, 'YLim', [min(freq0(1:1601,4))*1.1, max(freq0(1:1601,4))*1.1]);
% ylabel(gca, 'Angle output (rad)')

%% amplitude and magnitude values - Eyeball vs Goertzel
% eyeball
amp_gain_angle = [0.032/4, 0.03/4, 0.03/4, 0.02/4, 0.018/4];
amp_gain_pos = [0.01/4, 0.015/4, 0.01/4, 0.005/4];

phase_angle = [pi/2, pi/2, pi/2+0.3, pi/2+0.5, pi/2+0.5];
phase_pos = [1.3, 1.4, 1.7, 1.7];

% Goerztel
Fs = 200;
L =  2101;
logs = zeros(L, 4, length(freqs));
logs(:,:,1) = freq0(1:L,:); 
logs(:,:,2) = freq1(1:L,:);
logs(:,:,3) = freq2(1:L,:);
logs(:,:,4) = freq3(1:L,:); 
logs(:,:,5) = freq4(1:L,:);

gainX_data = zeros(1,length(freqs));
gainTheta_data = zeros(1,length(freqs));
phaseTheta_data = zeros(1,length(freqs));
phaseX_data = zeros(1,length(freqs));

for i = 1:length(freqs)
    log_temp = logs(:,:,i);
    fid = round( freqs(i)*L/(2*pi*Fs))+1;
    U = goertzel(log_temp(1:2101,1)', fid);
    pos = goertzel(log_temp(1:2101,3)', fid);
    angleTemp = goertzel(log_temp(1:2101,4)', fid);
    
    gainX_data(i) = abs(pos/U);
    gainTheta_data(i) = abs(angleTemp/U);
    phaseX_data(i) = angle(pos/U);
    phaseTheta_data(i) = angle(angleTemp/U);
end

fin_gainX = gainX_data; %amp_gain_pos
fin_gainTheta = gainTheta_data; %amp_gain_angle
fin_phaseX = phaseX_data; %phase_pos
fin_phaseTheta = phaseTheta_data; %phase_angle
beginPos = 1;% 2

%% plots for bode comparison
[magX,phaseX] = bode(cltf_X, {freqs(1)-1, freqs(end)+1});
magX = squeeze(magX); phaseX = squeeze(phaseX);
wradX = linspace(freqs(1)-1, freqs(end)+1, length(magX));

figure(1)
subplot(2,1,1)
semilogx(wradX, magX)
hold on
plot(freqs(beginPos:end),fin_gainX, 'r*')
hold off
grid
subplot(2,1,2)
semilogx(wradX, phaseX)
hold on
plot(freqs(beginPos:end), 180/pi*fin_phaseX, 'r*')
hold off
grid

[magTheta,phaseTheta] = bode(cltf_theta, {freqs(1)-1, freqs(end)+1});
magTheta = squeeze(magTheta); phaseTheta = squeeze(phaseTheta);
wradTheta = linspace(freqs(1)-1, freqs(end)+1, length(magTheta));

figure(2)
subplot(2,1,1)
semilogx(wradTheta, magTheta)
hold on
plot(freqs,fin_gainTheta, 'r*')
hold off
grid
subplot(2,1,2)
semilogx(wradTheta, phaseTheta)
hold on
plot(freqs, 180/pi*fin_phaseTheta, 'r*')
hold off
grid

figure(3)
bode(cltf_X)

figure(4)
bode(cltf_theta)

%% 