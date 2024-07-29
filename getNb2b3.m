function [n_1_lambda_array, GVD_array, TOD_array] = getNGVDTODfromLorentzLorenz(lambda_array, rho_medium, T_medium)

c_velocity = 3e8;
% rho_medium = 1000; %1000 kg/m^3
% lambda_medium = 532;
% T_medium = 293.15;
%%Lorentz-Lorenz Scattering
a_0 = 0.243905091;
a_1 = 9.53518094e-3;
a_2 = -3.64358110e-3;
a_3 = 2.65666426e-4;
a_4 = 1.59189325e-3;
a_5 = 2.45733798e-3;
a_6 = 0.897478251;
a_7 = -1.63066183e-2;

lambda_UV = 0.2292020;
lambda_IR = 5.432937;
rho_0 = 1000;
lambda_0 = 589e-9;  %589nm
T_0 = 273.15; %0 摄氏度
% rho = rho_medium / rho_0;
% T = T_medium / T_0;

syms n lambda T rho

    
    eq = a_0 + a_1 * (rho/ rho_0) + a_2 * (T / T_0) + a_3 * (lambda / lambda_0)^2 * (T / T_0) + ...
        a_4 / (lambda / lambda_0)^2 + a_5 / ((lambda / lambda_0)^2 - lambda_UV^2) + ...
        a_6 / ((lambda / lambda_0)^2 - lambda_IR^2) + a_7 * (rho/ rho_0)^2 - ...
        (n^2 - 1) / (n^2 + 2) / (rho/ rho_0) == 0;
    n = solve(eq, n);
%% 解参数方程得到两个n，留下大于0的n,n_1大于0
% n_1 = n(1);
% n_2 = n(2);
% lambda_medium = 532e-9;
% T_medium = 293.15;
% rho_medium = 1000; %1000 kg/m^3
% n_1 = subs(n_1, lambda, lambda_medium);
% n_2 = subs(n_2, lambda, lambda_medium);
% n_1 = subs(n_1, T, T_medium);
% n_2 = subs(n_2, T, T_medium);
% n_1 = subs(n_1, rho, rho_medium);
% n_2 = subs(n_2, rho, rho_medium);
% n_1 = double(n_1);
% n_2 = double(n_2);
% T_medium = T_0 + 24;
% rho_medium = 1000; %1000 kg/m^3
% lambda_532 = 532e-9;
n_parameter = n(1);  %%折射率的参数表示
n = n_parameter;
% n_diff_2 = diff(n, lambda, 2);

n_1_lambda = subs(n, T, T_medium);
n_1_lambda = subs(n_1_lambda, rho, rho_medium);
n_1_lambda_array = [];
% lambda_array = 400e-9 : 1e-9 : 1100e-9;
for i = 1 : 1 : length(lambda_array)
    n_1_lambda_array(i) = double(subs(n_1_lambda, lambda, lambda_array(i)));
end
% figure
% plot(lambda_array, n_1_lambda_array)
% title('reflective index')

n_diff_2 = diff(n_1_lambda, 2);
n_diff_3 = diff(n_1_lambda, 3);
% n_diff_2_value = double(subs(n_diff_2, lambda, lambda_532));
% n_diff_3_value = double(subs(n_diff_3, lambda, lambda_532));
GVD = lambda^3 / (2 * pi * c_velocity^2) * n_diff_2;
TOD = lambda^4 / (4 * pi^2 * c_velocity^3) * (3 * n_diff_2+ lambda * n_diff_3);
% GVD = GVD * 1e27;
% TOD = TOD * 1e42;
for i = 1 : 1 : length(lambda_array)
    GVD_array(i) = double(subs(GVD, lambda, lambda_array(i)));
    TOD_array(i) = double(subs(TOD, lambda, lambda_array(i)));
end

% 
% 
% GVD_array = [];
% TOD_array = [];
% 
% for i = 1 : 1 : length(lambda_array)
%     GVD_array(i) = double(subs(GVD, lambda, lambda_array(i))) ;
%     TOD_array(i) = double(subs(TOD, lambda, lambda_array(i)));
% end
% figure
% plot(lambda_array * 1e9, GVD_array)
% title('GVD')
% figure
% plot(lambda_array * 1e9, TOD_array)
% title('TOD')
% % 
% index_520 = find(lambda_array > 519e-9 & lambda_array < 521e-9);
% fprintf('5520nm, n = %f, GVD = %f, TOD = %f\n', n_1_lambda_array(index_520), GVD_array(index_520) * 1e27, TOD_array(index_520)* 1e42);
