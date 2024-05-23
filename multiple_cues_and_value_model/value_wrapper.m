num_meta_tr = 10;
value_list = [.25 .5 .75 1 1.25 1.5];
ff_values = zeros(length(value_list),num_meta_tr);
for valuei = 1:length(value_list)
    value = value_list(valuei);
    for tempi = 1:num_meta_tr
        run('Main_vectorized_simple_dopa.m')
        ff_values(valuei,tempi) = ff_vect(num_trials);
    end
end
%%
% figure;
% scatter(zeros(1,num_meta_tr),ff_values);

