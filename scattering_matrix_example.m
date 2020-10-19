%% program to create the scattering matrix
%
% the code works to create the fourier expansion and then to calculate the
% scattering matrix. This is because the calculation of the fourier
% expansion takes the longest time and requires little storage while the
% scattering matrix is a simple calculation but it large to store. 

radius = [5,10,25,50]; %  floe radius
floe_thickness = [0.5,1]; % floe thickness
period = [5:5:20]; % period

Fourier_cosine_D = Fourier_cosine_D_make(radius, floe_thickness,period);

%% test using the identies
% i = 34; % a random value
% b_n = Fourier_cosine_D(i).D_fourier;
% K = (2*pi/Fourier_cosine_D(i).period)^2/9.81;%wavenumber
% 
% abs(1 - sqrt(2*pi*K)*b_n(1))
% for n = 1:length(b_n)-1
%     abs(1 - (-1)^(n)*sqrt(pi*K/2)*b_n(n+1))
% end

%%

% we need to find the index for a given case
% e.g. T = 10s, floe_thickness = 1; radius = 25;

i = find([Fourier_cosine_D.period]== 10 ...
    & [Fourier_cosine_D.floe_thickness]== 1 & [Fourier_cosine_D.radius]== 25);

number_of_angles = 5;concentration = 0.3;time_step = 600;
% angle is evenly spaced between 0 and 2*pi. 

[S,exp_S] = scattering_matrix(Fourier_cosine_D(i),number_of_angles,concentration,time_step);

% the solution at going from time t to time t+time_step is given by
% multiplying by S