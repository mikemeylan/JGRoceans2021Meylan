function Fourier_cosine_D = Fourier_cosine_D_make(radius, floe_thickness,period)
%
% program to create the data file of Fourier_cosine expansion of D

if nargin == 0
radius = [5,10,25,50]; %  floe radius
floe_thickness = [0.5,1]; % floe thickness
period = [5:1:20]; % period
end

% these two variables we might want to change
E = 6e9; % Youngs modulus,
N=10; % number of Fourier modes

% create the structure which will store the data
Total_number = length(radius)*length(floe_thickness)*length(period);
Fourier_cosine_D = repmat( struct('radius',[], 'floe_thickness',[], 'period',[], 'D_fourier',[]), 1, Total_number);

%%
count = 1;
for count_radius = 1:length(radius)
    for count_floe_thickness = 1:length(floe_thickness)
        for count_period = 1:1:length(period)
            Fourier_cosine_D(count).radius = radius(count_radius);
            Fourier_cosine_D(count).floe_thickness = floe_thickness(count_floe_thickness);
            Fourier_cosine_D(count).period = period(count_period);
            g = 9.81;
            wave_number = (2*pi/period(count_period))^2/g;
            H = 2*pi/wave_number; % water depth set to be infinite while not too big
            omega = 2*pi/period(count_period);
            number_of_roots = 10; % used in the eigenfunction expansion.
            
            rho_w = 1000; % density of the water,
            rho_i = 900; % density of the ice
            nu = 0.3; % Poissons ratios
        
            D_coefficents = elastic_plate_D_coefficients(radius(count_radius),...,
                floe_thickness(count_floe_thickness),H,omega,E,rho_w,rho_i,nu,number_of_roots,N);
            %     D_2 = D_coefficents(1);
            %     for n = 1:N
            %         D_2 = D_2 + D_coefficents(n+1) * cos(n*theta);
            %     end
            
            Fourier_cosine_D(count).D_fourier = D_coefficents;
            count = count + 1;
        end
    end
end



