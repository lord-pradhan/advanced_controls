close all
sample_time = 0.005;
k_p = 0.01;

time_ms = time_ms - time_ms(1);
figure(1)
plot(time_ms/1e3, omega)
hold on
plot(out.stateOut)
legend('Hardware', 'Simulink')

figure(2)
plot(time_ms/1e3, control)
hold on
plot(out.control)
legend('Hardware', 'Simulink')