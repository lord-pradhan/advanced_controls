clear

sample_time = 0.1;
k_p = .17; K = 20.133; tau = 0.0985;
s = tf('s');
sysx= pade(exp(-sample_time*s/2),1);
G = k_p * K/(tau *s);
G_comb = G*sysx;
T = G_comb/(1+G_comb);
P = pole(T);