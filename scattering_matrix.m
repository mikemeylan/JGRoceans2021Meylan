function [S,exp_S_delta_t,S_alpha] = scattering_matrix(Fourier_cosine_D,number_of_angles,concentration,time_step)
% creates the scattering matrix
%
% Fourier_cosine_D is the coefficents of D(theta) as a fourier cosine
% series
% number of angles is the number of (different) angles.
% concentration is the concentration of floes between zero and one
% time_step is the time step in seconds

% we find D
D_coefficents = Fourier_cosine_D.D_fourier;
r = Fourier_cosine_D.radius;
period = Fourier_cosine_D.period;

c_g =  9.81*period/(4*pi);% the group_speed

theta = linspace(0,2*pi-2*pi/number_of_angles,number_of_angles);
D = D_coefficents(1);
for n = 2:length(D_coefficents)
    D = D + D_coefficents(n) * cos((n-1)*theta);
end

% we create the first row of S

S_row = c_g*2*pi/number_of_angles/(pi*r^2)*abs(D).^2;
S_row(1) = - sum(S_row(2:end));
S  = zeros(number_of_angles,number_of_angles);

for count = 1:number_of_angles
    S(:,count) = circshift(S_row,[0,count-1]);
end

exp_S_delta_t = expm(concentration*time_step*S);

S_alpha = -eye(size(S))*sum(c_g*2*pi/number_of_angles/(pi*r^2)*abs(D).^2);

%% you can calculate this matrix exponential really easily directly from S_row
%% This uses properties of circulant matrices https://en.wikipedia.org/wiki/Circulant_matrix
%
% uncomment below
%
% N_roots_unity = exp(2i*pi*[0:number_angles-1]/number_angles).';
% P = [];
% for count = 0:number_angles-1
%     P = [P,1/sqrt(number_angles)*N_roots_unity.^count]
% end
% 
% D = S_row(1); 
% 
% for count = 1:number_angles-1
%     D = D + S_row(number_angles-count+1)*N_roots_unity.^count;
% end
% 
% S_out2 = real(P*diag(exp(concentration*time_step*D)));