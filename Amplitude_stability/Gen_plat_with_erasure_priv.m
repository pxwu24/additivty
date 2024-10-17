a_plat = 3;
b_plat = 3;
e_plat = 3;

p = 0.5;
a_eras = 4;
b_eras = 5;
e_eras = 5;
V_eras = erasure(a_eras, p, 'isom');

dim = [a_plat * a_eras, b_plat * b_eras, e_plat * e_eras, 2 * a_plat * a_eras];
k = 2 * a_plat * a_eras;

num_points = 1;
sim_res = nan(num_points * 2, 3);
sim_res(:, 1) = .5 + .5 * rand(num_points * 2, 1);
sim_res(:, 2) = .5 * rand(num_points * 2, 1);

% Loop over possible values of s and t
for j = 1:1

    s = sim_res(j, 1);
    t = sim_res(j, 2);

    s = 0.3;
    t = 0.4;
    
    if s + t < 1

        V_plat = Generalized_vs_channel(s, t, 'isom'); % isomotry of generalized Platypus
        V = kron(V_plat, V_eras);
        V_2 = V;
        for k = 1:size(V_2,2)
            V_2(:, k) = syspermute(V_2(:, k), [1,3,2,4], [b_plat, e_plat, b_eras, e_eras]);
        end

        temp_private = optimize_private_information(V_2, dim, k);
        disp(['The parameters of generalized Platypus channel are: ', 's= ', num2str(s),', t= ', num2str(t)])
        disp(['private information of generalized Platypus channel = ', num2str(temp_private)])
        sim_res(j, 3) = temp_private;

    end
end