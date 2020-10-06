clc
clear all

syms z T
tau = 0.0985;
K = 20.133;
k_p = 2.12;
k_i = k_p/tau;
C = k_p + k_i*T/(z-1);
G = K*(1-exp(-T/tau))/(z-exp(-T/tau));
T = G*C/(1 + G*C);
[N,D] = numden(T);
eqn = D == 0;
X = solve(eqn);
pretty(vpa(X,2))
eqn2 =  X(2,1)+1;
eq2f = matlabFunction(eqn2);
x0 = 0.005;
fsolve(eq2f, 0.005)
