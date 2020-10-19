%% program to create a table of attenuation from the 3D scattering

clear all

concentration = 1;
time_step = 600;
number_of_angles = 24;% angle is evenly spaced between 0 and 2*pi.
radius = [10,20,50,100]; %  floe radius
floe_thickness = [0.5,1,1.5,2]; % floe thickness

% radius = [25]; %  floe radius
% floe_length = 10*randn(1,10);
% floe_length = 2*radius + floe_length - mean(floe_length);
floe_thickness = 0.1:0.1:2; % floe thickness

frequency = 0.03390909*((1.1).^[0:20]);
period = 1./frequency; % period

period = linspace(5,20,10);

number_of_angles = 100;

c_g =  9.81*period.'/(4*pi);
%%
%
% S = [];


%
output = zeros(length(period),length(floe_thickness));
output_2 = output;output_3 = output;
c_store = zeros(2,length(floe_thickness),length(radius));

for count_radius = 1:length(radius)
    
    for count_thickness = 1:length(floe_thickness)
        count_thickness
        
        
        Fourier_cosine_D = Fourier_cosine_D_make(radius(count_radius), floe_thickness(count_thickness),period);
        
        
%         R_and_T = R_and_T_make(radius(count_radius), floe_thickness(count_thickness),period);
        for i = 1:length(Fourier_cosine_D)
            
%             for count2 = 1:length(floe_length)
%                 transmission = transmission_make(floe_length(count2),...
%                     floe_thickness(count_thickness),period(i));
%                 T_store(count2) = transmission.T;
%             end
            [S,~] = scattering_matrix(Fourier_cosine_D(i),number_of_angles,concentration,time_step);
            c_g =  9.81*Fourier_cosine_D(i).period/(4*pi);
            output(i,count_thickness) = -S(1,1)/c_g;
%             T = R_and_T(i).R_and_T(2);
%             output_2(i,count_thickness) = -log(abs(T^2))/(2*radius);
%             output_3(i,count_thickness) = -log(mean(abs(T_store.^2)))/(2*radius);
        end
        
        
        
        %%
        c = [log(period).',ones(size(period)).']\(log(output(:,count_thickness)));
        
        % plot(log(period),log(output),log(period),c(1)*log(period) + c(2))
        % semilogy(period,output,period,exp(c(2))*period.^c(1),period,output_2)
        % hold on
        % semilogy(period,-log(mean(abs(T_store)).^2)/mean(floe_length))
        
        c_store(:,count_thickness,count_radius) = c;
        
    end
    
end

% a = 2.12e-3; b = 4.59e-2;
%
% alpha_2 = a*period.^(-2) + b*period.^(-4);
%
% semilogy(period,alpha_2)
%
% hold off

%% horvat
i=5;
log_alpha = -0.3203 + 2.058*floe_thickness - 0.9375*period(i) ...
    - 0.4269*floe_thickness.^2 + 0.1566*floe_thickness*period(i) + ....
    0.0006*period(i).^2 - log(2*radius);


plot(floe_thickness,log(output(i,:)),floe_thickness,log(output_2(i,:))...
    ,floe_thickness,log_alpha,floe_thickness,log(output_3(i,:)))


%%

%% horvat
i=10;
log_alpha = -0.3203 + 2.058*floe_thickness(i) - 0.9375*period ...
    - 0.4269*floe_thickness(i).^2 + 0.1566*floe_thickness(i)*period + ....
    0.0006*period.^2 - log(2*radius);


plot(period,log(output(:,i)),period,log(output_2(:,i))...
    ,period,log_alpha,period,log(output_3(:,i)))
