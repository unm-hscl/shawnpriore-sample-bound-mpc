%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clean environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all;
cvx_clear;

quiet = 0;

problem_setup;

solve_our_method;
solve_cantelli;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%
colors = [224,   0,   0; % red
           30, 144  255; % dark blue
            0, 170,  85; % green
          118,   0, 168; % purple
           46,  52,  59; % grey
          236, 176,  31;  % yellow
           76, 189, 237; % light blue
          161,  19,  46  % dark red
           ] ./ 255;
       
plot_symbols = ['o', 'd', '^', 'h', 'v', '>', 'p', 's'];

Faces = [1,2,3,4];

fig = figure();

subplot(8,1,1)
hold on

plot(nan, nan, '--k', 'Marker', 'none');
plot(nan, nan, '-k', 'Marker', 'none');

patch('Faces',Faces,'Vertices', Polyhedron([-eye(2); eye(2)], ones(4,1)).V,...
    'FaceColor',  'white', ...
    'EdgeColor', 'black', ...
    'FaceAlpha', 0); 

plot(nan, nan, 'Marker', 's', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineStyle', 'none');


for i = 1:3
   plot(nan, nan, 'Color', colors(i, :), 'Marker', plot_symbols(i));
end

plots=get(gca, 'Children');

legend([plots(3), plots(2), plots(1), plots(5), plots(4), plots(6), plots(7)], ...
     {'SAT 1', 'SAT 2', 'SAT 3', 'Target Set', 'Initial Location', 'Proposed', 'MPC With Cantelli' },...
    'Orientation','horizontal', ...
    'Location', 'south', ...
    'NumColumns', 7, ...
    'interpreter', 'latex');

axis([0 0.1 0 0.1]);
axis off
hold off


subplot(8,1,[2:4]);
hold on

patch('Faces',Faces,'Vertices', Polyhedron(target_set_a.A([1;2;7;8],1:2), target_set_a.b([1;2;7;8])).V,...
    'FaceColor', colors(1,:), ...
    'FaceAlpha', 0.1); 
patch('Faces',Faces,'Vertices', Polyhedron(target_set_b.A([1;2;7;8],1:2), target_set_b.b([1;2;7;8])).V,...
    'FaceColor', colors(2,:), ...
    'FaceAlpha', 0.1); 
patch('Faces',Faces,'Vertices', Polyhedron(target_set_c.A([1;2;7;8],1:2), target_set_c.b([1;2;7;8])).V,...
    'FaceColor', colors(3,:), ...
    'FaceAlpha', 0.1); 

plot( [x_0_a(1); Ex_A(1:6:end)], [x_0_a(2); Ex_A(2:6:end)], '-', 'Color', colors(1,:), 'Marker', plot_symbols(1));
plot( [x_0_b(1); Ex_B(1:6:end)], [x_0_b(2); Ex_B(2:6:end)], '-', 'Color', colors(2,:), 'Marker', plot_symbols(2));
plot( [x_0_c(1); Ex_C(1:6:end)], [x_0_c(2); Ex_C(2:6:end)], '-', 'Color', colors(3,:), 'Marker', plot_symbols(3));

plot( [x_0_a(1); Ex_A_r(1:6:end)], [x_0_a(2); Ex_A_r(2:6:end)], '--', 'Color', colors(1,:), 'Marker', plot_symbols(1));
plot( [x_0_b(1); Ex_B_r(1:6:end)], [x_0_b(2); Ex_B_r(2:6:end)], '--', 'Color', colors(2,:), 'Marker', plot_symbols(2));
plot( [x_0_c(1); Ex_C_r(1:6:end)], [x_0_c(2); Ex_C_r(2:6:end)], '--', 'Color', colors(3,:), 'Marker', plot_symbols(3));



plot( x_0_a(1), x_0_a(2), 'k', 'MarkerFaceColor', 'k', 'Marker', plot_symbols(1));
plot( x_0_b(1), x_0_b(2), 'k', 'MarkerFaceColor', 'k', 'Marker', plot_symbols(2));
plot( x_0_c(1), x_0_c(2), 'k', 'MarkerFaceColor', 'k', 'Marker', plot_symbols(3));

xlabel('x');
ylabel('y');
hold off

subplot(8,1,[6:8]);
hold on

patch('Faces',Faces,'Vertices', Polyhedron(target_set_a.A([1;3;7;9],[1;3]), target_set_a.b([1;3;7;9])).V,...
    'FaceColor', colors(1,:), ...
    'FaceAlpha', 0.1); 
patch('Faces',Faces,'Vertices', Polyhedron(target_set_b.A([1;3;7;9],[1;3]), target_set_b.b([1;3;7;9])).V,...
    'FaceColor', colors(2,:), ...
    'FaceAlpha', 0.1); 
patch('Faces',Faces,'Vertices', Polyhedron(target_set_c.A([1;3;7;9],[1;3]), target_set_c.b([1;3;7;9])).V,...
    'FaceColor', colors(3,:), ...
    'FaceAlpha', 0.1); 

plot( [x_0_a(1); Ex_A(1:6:end)], [x_0_a(3); Ex_A(3:6:end)], '-', 'Color', colors(1,:), 'Marker', plot_symbols(1));
plot( [x_0_b(1); Ex_B(1:6:end)], [x_0_b(3); Ex_B(3:6:end)], '-', 'Color', colors(2,:), 'Marker', plot_symbols(2));
plot( [x_0_c(1); Ex_C(1:6:end)], [x_0_c(3); Ex_C(3:6:end)], '-', 'Color', colors(3,:), 'Marker', plot_symbols(3));

plot( [x_0_a(1); Ex_A_r(1:6:end)], [x_0_a(3); Ex_A_r(3:6:end)], '--', 'Color', colors(1,:), 'Marker', plot_symbols(1));
plot( [x_0_b(1); Ex_B_r(1:6:end)], [x_0_b(3); Ex_B_r(3:6:end)], '--', 'Color', colors(2,:), 'Marker', plot_symbols(2));
plot( [x_0_c(1); Ex_C_r(1:6:end)], [x_0_c(3); Ex_C_r(3:6:end)], '--', 'Color', colors(3,:), 'Marker', plot_symbols(3));


plot( x_0_a(1), x_0_a(3), 'k', 'MarkerFaceColor', 'k', 'Marker', plot_symbols(1));
plot( x_0_b(1), x_0_b(3), 'k', 'MarkerFaceColor', 'k', 'Marker', plot_symbols(2));
plot( x_0_c(1), x_0_c(3), 'k', 'MarkerFaceColor', 'k', 'Marker', plot_symbols(3));

xlabel('x');
ylabel('z');
hold off