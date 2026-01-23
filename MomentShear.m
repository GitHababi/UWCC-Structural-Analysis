close all;
clc;

% Run BuoyantLoad.m before this

%% LOADS
%Point Loads
paddlers = P
self_weight = selfWeight;
buoyant_force = buoyantLoad;

%% SHEAR AND MOMENT CALCS

% init vars
x_common = linspace(min(self_weight(:,1)),max(self_weight(:,1)),1000);
load = buoyant_force(:,2) + self_weight(:,2);

load_interpolated = interp1(buoyant_force(:,1), load,x_common,'linear');
V =  cumtrapz(x_common, load_interpolated);

% adjust shear for point loads
for i = 1:size(paddlers,1)
    for j = 1:1000
        if x_common(j) >= paddlers(i,1)
            V(j) = V(j) + paddlers(i,2);
        end
    end
end

M = cumtrapz(x_common, V);


%% PLOT

hold on
plot(x_common,V,'Color',[0/255,100/255,100/255],'LineWidth',2,'DisplayName','Shear');
xlim([x_common(1) x_common(end)]);
title('Shear Diagram');
xlabel('x-position or length of canoe (in.)');
ylabel('Shear (lbs.)');

figure;

plot(x_common,M,'Color', [235/255 157/255 80/255], 'LineWidth', 2, 'DisplayName','Moment');
xlim([x_common(1) x_common(end)]);
title('Moment Diagram');
xlabel('x-position or length of canoe (in.)');
ylabel('Bending Moment (lbs.*in.)');
hold off

[V_m, Ivm] = max(V);
[V_mn, Ivmn] = min(V);
[M_m, Mvm] = max(M);
[M_mn, Mvmn] = min(M);
disp("V_max pos: " + V_m + "(lbs) @ " + x_common(Ivm))
disp("V_max neg: " + V_mn + "(lbs) @ " + x_common(Ivmn))
disp("M_max pos: " + M_m + "(lbs*in) @ " + x_common(Mvm))
disp("M_max neg: " + M_mn + "(lbs*in) @ " + x_common(Mvmn))