% Set up timings
n_procs = [1,2,4,8,16,32,64];
t_par_smooth = [.41,.41,.43,.44,.45,.66,.86];

t_ser_smooth = zeros(1,numel(n_procs));

for i = 1:1:numel(n_procs)
    t_ser_smooth(i) = t_par_smooth(1)*n_procs(i);
end

speedup = t_ser_smooth./t_par_smooth;
rel_dif = zeros(1,numel(speedup));
for i = 1:1:numel(speedup)
    rel_dif(i) = (speedup(i) - n_procs(i))/n_procs(i);
end
rel_dif

clf;
figure(1);
hold on;
grid on;
plot(n_procs,speedup,'ro-');
plot(n_procs,n_procs);
    legend('smoothing speedup', 'ideal speed up');
    xlabel('procs');
    ylabel('speedup');
    title('Weak scaling speedup of smoothing kernel'); 