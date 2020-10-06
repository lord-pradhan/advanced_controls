clear 
clc
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

C = [0,0,1,0];
%     0,0,0,1]; 
% C = eye(4);
D = zeros(size(C,1), size(B,2));

[n,d] = ss2tf(A,B,C,D);
%% pole placement
% k = place(A, B, [-0.01, -0.4, -35, -80])

%% lqr
Q = diag([0.001, 0.1, 5000, 1]);
R = 0.2;

[k,S,CLP] = lqr(A,B,Q,R)

%gains that work - {1.03, -113, -6.2}, {1.03, -154.6, -3}