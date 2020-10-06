clear

sample_time = 0.005; tau = 0.0985;

k_p = 2.5567; k_i = 5.444; k_d = 0.2217;
N = 218;
% 
a_0 = -1.5; a_1 = 0.5;
b_0 = k_p + 100*k_d; b_1 = -1.5*k_p + k_i/200 - 200*k_d; 
b_2 = 100*k_d - k_i/400 + k_p/2;

% syms z k_p k_i k_d
% T = k_p + k_i*sample_time/(z-1) + k_d*(N)/(1 + N*sample_time/(z-1) );
% [N, D] = numden(T);
% pretty(vpa(expand(N),3))
% pretty(vpa(expand(D),3))
% pretty(vpa(simplify(T),3))

% k_p = 0.05; k_i = k_p/tau;
% b_0 = k_p; b_1 = k_i*sample_time - k_p; a_1 = -1;