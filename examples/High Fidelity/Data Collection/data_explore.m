data = readmatrix('a_dist.dat');

fig = figure();

subplot(2,3,1);
hold on
histogram(data(1,:),20,'Normalization','pdf');
xlabel('$x$', 'interpreter', 'latex');
ylabel('$f(x)$', 'interpreter', 'latex');
hold off

subplot(2,3,2);
hold on
histogram(data(2,:),20,'Normalization','pdf');
xlabel('$y$', 'interpreter', 'latex');
ylabel('$f(y)$', 'interpreter', 'latex');
hold off

subplot(2,3,3);
hold on
histogram(data(3,:),20,'Normalization','pdf');
xlabel('$z$', 'interpreter', 'latex');
ylabel('$f(z)$', 'interpreter', 'latex');
hold off

subplot(2,3,4);
hold on
histogram(data(4,:),20,'Normalization','pdf');
xlabel('$\dot x$', 'interpreter', 'latex');
ylabel('$f(\dot x)$', 'interpreter', 'latex');
hold off

subplot(2,3,5);
hold on
histogram(data(5,:),20,'Normalization','pdf');
xlabel('$\dot y$', 'interpreter', 'latex');
ylabel('$f(\dot y)$', 'interpreter', 'latex');
hold off

subplot(2,3,6);
hold on
histogram(data(6,:),20,'Normalization','pdf');
xlabel('$\dot z$', 'interpreter', 'latex');
ylabel('$f(\dot z)$', 'interpreter', 'latex');
hold off