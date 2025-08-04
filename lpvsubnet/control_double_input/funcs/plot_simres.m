function plot_simres(res)
figure()
subplot(411)
plot(res.q4.Time, res.q4.Data); hold on
plot(res.q4ref.Time, res.q4ref.Data)
grid on
legend('$q_4$', '$q_{4ref}$', 'interpreter', 'latex')

subplot(412)
plot(res.q4.Time, res.q4.Data); hold on
plot(res.hq4.Time, res.hq4.Data)
grid on
legend('$q_4$', '$\hat q_4$', 'interpreter', 'latex')

subplot(413)
plot(res.i2.Time, res.i2.Data)
grid on
ylabel('$i_2$', 'interpreter', 'latex')

subplot(414)
plot(res.q4dot.Time, res.q4dot.Data); hold on
plot(res.hq4dot.Time, res.hq4dot.Data)
grid on
legend('$\dot q_4$', '$\dot{\hat q_4}$', 'interpreter', 'latex')
end
