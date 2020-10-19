%% program to create the scattering matrix table
%
% the code works to create the fourier expansion and then to calculate the
% scattering matrix. This is because the calculation of the fourier
% expansion takes the longest time and requires little storage while the
% scattering matrix is a simple calculation but it is large to store.
clear all
clc

concentration = 0.6;
time_step = 600;
number_of_angles = 24;% angle is evenly spaced between 0 and 2*pi.

% number_of_angles = 8;

radius = [10]; %  floe radius
floe_thickness = [1.5]; % floe thickness
frequency = 0.03390909*((1.1).^[0:32]);
period = 1./frequency; % period
theta = linspace(0,2*pi-2*pi/number_of_angles,number_of_angles);

period_count = 12; % around 10 s 

Fourier_cosine_D = Fourier_cosine_D_make(radius, floe_thickness,period(period_count));


%%

    [S,exp_S_delta_t] = scattering_matrix(Fourier_cosine_D(1),number_of_angles,concentration,time_step);
  
c_g =  9.81*period.'/(4*pi);

beta = exp_S_delta_t(1:number_of_angles:end,1);

alpha = ((1-beta)/time_step)./c_g(period_count);

D = diag(c_g(period_count)*cos(theta));

% D = diag(cos(theta));

% S = S/abs(S(1,1))

[V,D1] = eig(S,D)

[y,i] = sort(real(diag(D1)))

lambda = [y(1:number_of_angles/2 -2);0]

V1 = [V(:,i(1:number_of_angles/2 -2)),ones(number_of_angles,1)]

% need the positive angles

positive_index = [1:number_of_angles/4,...
    number_of_angles-number_of_angles/4+2:number_of_angles]

input = zeros(size(positive_index.'));input(1) = 1;

V2 = V1(positive_index,:);

c = V2\input

x = linspace(0,1000000,100);

output = 0;

for count = 1:length(c)
    output = output + (c(count)*V1(:,count))*exp(x*lambda(count));
end

lambda = [y(1:number_of_angles/2 -2);max(y(1:number_of_angles/2 -2))]

output2 = 0;

for count = 1:length(c)
    output2 = output2 + (c(count)*V1(:,count))*exp(x*lambda(count));
end

a = 2.12e-3; b = 4.59e-2;

alpha_2 = a*period(period_count).^(-2) + b*period(period_count).^(-4);

plot(x/1000,output(1,:),x/1000,exp(x*S(1,1)/c_g(period_count)),x/1000,output2(1,:),...
    x/1000,exp(-alpha_2*x))

