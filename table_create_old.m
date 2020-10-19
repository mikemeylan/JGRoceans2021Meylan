%% program to create the scattering matrix table
%
% the code works to create the fourier expansion and then to calculate the
% scattering matrix. This is because the calculation of the fourier
% expansion takes the longest time and requires little storage while the
% scattering matrix is a simple calculation but it is large to store.
clear all
warning('off','all')

concentration = 0.6;
time_step = 600;
number_of_angles = 24;% angle is evenly spaced between 0 and 2*pi.
radius = [5, 10] %, 15, 20, 25, 30, 35, 40, 45, 50 ]; %  floe radius
floe_thickness = [0.5, 1.0] % , 1.5, 2.0, 2.5, 3.0]; % floe thickness
frequency = 0.03390909*((1.1).^[0:32]);
period = 1./frequency; % period

Fourier_cosine_D = Fourier_cosine_D_make(radius, floe_thickness,period);


%%

S = [];
S2 = [];
S3 = [];
S4 = [];
S_matrix = [];

for i = 1:length(period)
    [S_orginal,exp_S_delta_t] = scattering_matrix(Fourier_cosine_D(i),number_of_angles,concentration,time_step);
    S = [S; exp_S_delta_t];
    S_matrix = [S_matrix; S_orginal];
    exp_S_delta_t2 = 0.5*diag(diag(exp_S_delta_t)) +0.5*exp_S_delta_t;
    % halves the terms off the diagonal.
    S2 = [S2;exp_S_delta_t2];
    
    exp_S_delta_t3 = 0.9*diag(diag(exp_S_delta_t)) +0.1*exp_S_delta_t;
    % 0.1 the terms off the diagonal.
    S3 = [S3;exp_S_delta_t3];
    
    exp_S_delta_t4 = 1*diag(diag(exp_S_delta_t)) +0*exp_S_delta_t;
    % zero the terms off the diagonal.
    S4 = [S4;exp_S_delta_t4];
end

%%

save('bash_test_data_600.dat','S','-ascii')

save('bash_test_data_600','S_matrix')

save('bash_test_data_600_energy_loss.dat','S2','-ascii')

save('bash_test_data_600_energy_loss3.dat','S3','-ascii')

save('bash_test_data_600_energy_loss4.dat','S4','-ascii')

c_g =  9.81*period.'/(4*pi);
beta = S(1:number_of_angles:end,1);

alpha = ((1-beta)/time_step)./c_g;

a = 2.12e-3; b = 4.59e-2;

alpha_2 = a*period.^(-2) + b*period.^(-4);


S5=[];
for i = 1:length(period)
    S5 = [S5; exp(-alpha_2(i)*time_step*c_g(i))*eye(number_of_angles,number_of_angles)];
    
end

S6=[];
for i = 1:length(period)
    S6 = [S6; 1/2*exp(-alpha_2(i)*time_step*c_g(i))*eye(number_of_angles,number_of_angles)];
    
end

beta2 = S5(1:number_of_angles:end,1);

alpha_3 = ((1-beta2)/time_step)./c_g;

semilogy(period,alpha,period,alpha_2,period,alpha_3)

save('bash_test_data_600_meylan_etal2014.dat','S5','-ascii')

save('bash_test_data_600_meylan_etal2014_halved.dat','S6','-ascii')
