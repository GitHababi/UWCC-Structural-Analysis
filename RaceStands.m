close all;
clc;

% Run BuoyantLoad.m before this

%% Change first column to choose race stand position.
stands = [
    -72 0;
     72 0
];


%% Calc Moment, Shear, Reaction Forces
self_weight = selfWeight;
x_common = linspace(min(self_weight(:,1)),max(self_weight(:,1)),1000);

load_interpolated = interp1(self_weight(:,1), self_weight(:,2),x_common,'linear');
V =  cumtrapz(x_common, load_interpolated);

dx = x_common(2) - x_common(1);
moment_stand1 = sum(load_interpolated .* dx .* (x_common - stands(1,1)));

span = stands(2,1) - stands(1,1);
stands(2,2) = moment_stand1 / span;
stands(1,2) = trapz(x_common, load_interpolated) - stands(2,2);

% adjust shear for stand loads
for i = 1:size(stands,1)
    for j = 1:1000
        if x_common(j) >= stands(i,1)
            V(j) = V(j) - stands(i,2);
        end
    end
end

M = cumtrapz(x_common, V);
%% PLOT Free Body Diagram

figure('Color','w'); hold on; grid on;
plot(xUnique, q_self, 'r--', 'LineWidth',1.6);
yline(0,'k-');

for i = 1:size(P,1)
    x0 = stands(i,1);
    Pmag = abs(stands(i,2));
    quiver(x0, 0, 0, Pmag/100, 0, 'MaxHeadSize',1, 'Color','k', 'LineWidth',1.5); % represent paddlers as arrows
end

ylim([min(q_self)*1.5, abs(min(q_self)) *1.5])

xlabel('x (in)');
ylabel('Load (lbf/in)');
title('Self Weight, and Race Stands');
legend('Self Weight (down)','Race Stands (up)','Location','best');

%% PLOT Moment And Shear

figure;
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