function transmission = transmission_make(floe_length, floe_thickness,period)
%
% program to create the data file of Fourier_cosine expansion of D
warning('off','MATLAB:nearlySingularMatrix' )

if nargin == 0
floe_length = [5,10,25,50]; %  floe radius
floe_thickness = [0.5,1]; % floe thickness
period = [5:1:20]; % period
end

% these two variables we might want to change
E = 6e9; % Youngs modulus,
N=10; % number of Fourier modes

% create the structure which will store the data
Total_number = length(floe_length)*length(floe_thickness)*length(period);
transmission = repmat( struct('floe_length',[], 'floe_thickness',[], 'period',[], 'T',[]), 1, Total_number);

%%
count = 1;
for count_length = 1:length(floe_length)
    for count_floe_thickness = 1:length(floe_thickness)
        for count_period = 1:1:length(period)
            transmission(count).floe_length = floe_length(count_length);
            transmission(count).floe_thickness = floe_thickness(count_floe_thickness);
            transmission(count).period = period(count_period);
            g = 9.81;
            wave_number = (2*pi/period(count_period))^2/g;
            H = 2*pi/wave_number; % water depth set to be infinite while not too big
            omega = 2*pi/period(count_period);
            number_of_roots = 10; % used in the eigenfunction expansion.
            
            rho_w = 1000; % density of the water,
            rho_i = 900; % density of the ice
            nu = 0.3; % Poissons ratios
        
            g = 9.81;

% code uses these non-dimensional variables
alpha = omega^2/g;
D = E*floe_thickness(count_floe_thickness)^3/(12*(1-nu^2));
beta = D/(g*rho_w);
gamma = rho_i*floe_thickness(count_floe_thickness)/rho_w;

             [a_s,~,a_a,~] = finite_plate_symmetry(alpha,beta,gamma,H,number_of_roots,floe_length(count_length),0,nu);
    
            transmission(count).T = 1/2*(a_s(1) - a_a(1));
            count = count + 1;
        end
    end
end



