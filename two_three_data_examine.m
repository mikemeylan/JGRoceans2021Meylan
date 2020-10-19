%%
clear all

load  two_three_data

radius_start = 50;
floe_thickness_start = 0.75;

%%
count_2 = 1;
for count  = 1:length(R_and_T)
    if (R_and_T(count).radius == radius_start) && ...
            (R_and_T(count).floe_thickness == floe_thickness_start)
        out(count_2,:) = R_and_T(count).R_and_T;
              out2(count_2,:) = Fourier_cosine_D(count).D_fourier;
        period(count_2) = R_and_T(count).period;
            floe_thickness(count_2,:) = R_and_T(count).floe_thickness;
        index(count_2) = count;
        count_2 = count_2 + 1;
    end
end
    
%%
loglog(period,-log(1-abs(out(:,1)).^2)/radius_start,...
  period,2e-3*period.^(-2),'g',period,-log(1-abs(out2(:,1)).^2)/radius_start,'r--')
% hold on
% loglog(period,-1e-4*period.^(-2))
% hold on



