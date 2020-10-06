clear all

sample_time = 0.005;
K = 20.133; tau= 0.0985;

s = tf('s');
G_prime = K/(s*(tau*s+1));
G_z = c2d(G_prime, sample_time);