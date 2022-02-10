figure

subplot(2,1,1)
plot(t1,Mby_1)
hold on
plot(t1,Mby_2)
hold on
plot(t1,Mby_3)
xlabel('t [s]')
ylabel('Mby [Nm]')
grid on

subplot(2,1,2)
plot(t1,smooth(Twr_FA,10))
xlabel('t [s]')
ylabel('thrust [N]')
grid on


figure

subplot(4,1,4)
plot(t1,Torque_f)
xlabel('t [s]')
ylabel('toruqe [Nm]')
ylim([-1, 5]);
grid on

subplot(4,1,2)
plot(t1,RPM_res_f)
xlabel('t [s]')
ylabel('rpm [1/min]')
ylim([0 600]);
grid on

subplot(4,1,3)
plot(pitch(:,1),sum(pitch(:,[2,4,6]),2))
xlabel('t [s]')
ylabel('pitch [deg]')
% ylim([0 600]);
grid on

subplot(4,1,1)
plot(t1,hw_v1)
xlabel('t [s]')
ylabel('pitch [deg]')
ylim([2 14]);
grid on
