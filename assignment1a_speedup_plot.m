% Set up timings
n_procs = [1,2,8,16];
t_par_count_x = [0.09, 0.08, 0.09, .10];
t_par_count_y = [0.08, 0.09, 0.09, .10];

t_ser_count_x = zeros(1,numel(n_procs));
t_ser_count_y = zeros(1,numel(n_procs));

for i = 1:1:numel(n_procs)
    t_ser_count_x(i) = t_par_count_x(1)*n_procs(i);
    t_ser_count_y(i) = t_par_count_y(1)*n_procs(i);
end

speedup_x = t_ser_count_x./t_par_count_x;
speedup_y = t_ser_count_y./t_par_count_y;

clf;
figure(1);
hold on;
grid on;
plot(n_procs,speedup_x);
plot(n_procs,speedup_y);
plot(n_procs,n_procs);
    legend('x speedup', 'y speed up', 'theoretical speed up');