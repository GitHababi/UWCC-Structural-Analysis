close all;
clc;

% Set paddler information (position, FOS, weight) in generatePaddlerLoads function.

%% Constants 
gamma = 62.4/1728;        % water density (lbf/in^3)
xref  = 0;
draft_in = 0.63*12;       % change to 0.82*12 for 4-paddler or .63*12 for male tandem (in)
dead_FOS = 1.5;
rho_struct = 68.5/1728;   % density of structural mix of concrete (lbf/in^3)
interp_rate = 1000;       % number of subdivided data points.

%% Load Paddler Cases
paddlerCases = generatePaddlerLoads();           
P = paddlerCases.case_4;  % Choose paddler case here.

%% Load Cross Section Information.
coordinates     = readtable('HullCoordinates.csv'); % input CSV
xCoordinates    = coordinates.("LongitudalPosition");
side            = string(coordinates.Side);
yCoordinates    = coordinates.Offset; 
zCoordinates    = coordinates.Height;

%% loop through stations and calculate CSA
xUnique = unique(xCoordinates); 
N = numel(xUnique);

csa = zeros(N,1); % CSA
sectionsOuter = cell(N,1); % cell is a matlab array that can contain diff data types
ymin = inf; ymax = -inf; zmin = inf; % bounds

for i = 1:N
    xi = xUnique(i);

    % OUTER
    mO = (xCoordinates==xi) & (side=="Outer");
    yO = yCoordinates(mO);  zO = zCoordinates(mO);

    % INNER
    mI = (xCoordinates==xi) & (side=="Inner");
    yI = yCoordinates(mI);  zI = zCoordinates(mI);

    [yO, idxO] = sort(yO,'descend'); zO = zO(idxO);

    yO_full = [yO; -flipud(yO)];
    zO_full = [zO;  flipud(zO)];

    yI_full = [yI; -flipud(yI)];
    zI_full = [zI;  flipud(zI)];

    y_loop = [ yO_full; flipud(yI_full) ];
    z_loop = [ zO_full; flipud(zI_full) ];

    dup = [false; abs(diff(y_loop)) + abs(diff(z_loop)) < 1e-9];
    y_loop(dup) = []; z_loop(dup) = [];

    csa(i) = abs(polyarea(y_loop, z_loop)); % polyarea!!

    sectionsOuter{i} = [yO_full, zO_full]; % store results to use later on
    ymin = min(ymin, min(yO_full)); ymax = max(ymax, max(yO_full)); 
    zmin = min(zmin, min(zO_full));                        
end

%% Buoyant Load with 0 trim and no paddlers (pure geometry)
theta = 0;                               
[~,k0] = min(abs(xUnique-0)); % closest x to midship 
z_keel0 = min(sectionsOuter{k0}(:,2)); % keel coordinate (lowest Z value)
z0 = z_keel0 + draft_in;   % z0 = keel + draft, first z thats above water               

padY = 10*max(1,(ymax-ymin)); padZ = 10*max(1,max(abs([ymin,ymax,zmin]))); % creates a rectangle representing water
A_sub = zeros(N,1);

% What this loop does: for each X position (cross section), overlays a
% rectangle where the top edge is the waterline (below the water)
for j = 1:N
    z_wl = z0 + (xUnique(j) - xref)*tan(theta);
    clipY = [ymin-padY, ymax+padY, ymax+padY, ymin-padY];
    clipZ = [zmin-padZ, zmin-padZ, z_wl,      z_wl     ];
    Pout  = polyshape(sectionsOuter{j}(:,1), sectionsOuter{j}(:,2), 'Simplify', true); % define shape
    Pwet  = intersect(Pout, polyshape(clipY, clipZ)); % finds intersection of rectangle and cross section, finds wetted section
    if ~isempty(Pwet.Vertices)
        A_sub(j) = area(Pwet); end % finds Area of wetted section
end
q_buoy = gamma * A_sub;       % standard buoyancy calc   
q_buoy_perCSA = q_buoy ./ csa;   % buoyant load per CSA

%% Preliminary printout
csas = table(xUnique, csa, A_sub, q_buoy, q_buoy_perCSA, ...
    'VariableNames', {'xStation_in','MaterialArea_in2','SubmergedArea_in2','BuoyantLoad_lbf_per_in','BuoyantLoad_perCSA'});
disp(csas);

figure; plot(xUnique, q_buoy, 'LineWidth',1.4); grid on;
xlabel('x (in)'); ylabel('Buoyant Load (lbf/in)'); title('Buoyant Load (θ=0 rad)');

%% Include self-weight + paddlers (equilibrium)
q_self = -rho_struct * csa;    % self weight load      
                      

[~,k0] = min(abs(xUnique-0));
z_keel0 = min(sectionsOuter{k0}(:,2));
v0 = [z_keel0 + 0.63*12, 0];  % inital waterline + trim, will update and solve for equillibrium 

ymin_all = min(cellfun(@(v)min(v(:,1)),sectionsOuter));
ymax_all = max(cellfun(@(v)max(v(:,1)),sectionsOuter));
zmin_all = min(cellfun(@(v)min(v(:,2)),sectionsOuter));
padY = 10*max(1,(ymax_all - ymin_all)); 
padZ = 10*max(1,max(abs([ymin_all,ymax_all,zmin_all])));

Ps = cellfun(@(v) polyshape(v(:,1),v(:,2),'Simplify',true), sectionsOuter, 'UniformOutput', false);
clipY = [ymin_all-padY, ymax_all+padY, ymax_all+padY, ymin_all-padY];

submerged = @(v) submerged_impl(v, Ps, xUnique, xref, gamma, clipY, zmin_all, padZ);          
resid     = @(v) resid_impl(v, submerged, xUnique, q_self, P, xref);                          

obj = @(vv) sum(resid(vv).^2);
v = fminsearch(obj, v0, optimset('Display','off')); % fminsearch searches for the v closest to 0
[z0, theta] = deal(v(1), v(2));
[A_sub, q_buoy] = submerged(v);

fprintf('Equilibrium: θ = %.4f rad (%.3f°), z₀ = %.3f in\n', theta, theta*180/pi, z0);

%% Final plot: buoyancy, self-weight, paddlers
figure('Color','w'); hold on; grid on;
plot(xUnique, q_buoy, 'b-', 'LineWidth',1.6);
plot(xUnique, q_self, 'r--', 'LineWidth',1.6);
yline(0,'k-');

for i = 1:size(P,1)
    x0 = P(i,1);
    Pmag = abs(P(i,2));
    quiver(x0, 0, 0, -Pmag/200, 0, 'MaxHeadSize',1, 'Color','k', 'LineWidth',1.5); % represent paddlers as arrows
end

xlabel('x (in)');
ylabel('Load (lbf/in)');
title('Buoyant Load, Self Weight, and Paddlers');
legend('Buoyancy (up)','Self Weight (down)','Paddlers (down)','Location','best');

%% Helper functions

% Calculates the buoyant load for a proposed trim 
function [A_sub,qb]=submerged_impl(v, Ps, xUnique, xref, gamma, clipY, zmin_all, padZ)
    Nloc = numel(Ps); A_sub = zeros(Nloc,1); % num different possible stations and areas
    for jj=1:Nloc
        z_wl = v(1) + (xUnique(jj)-xref)*tan(v(2));
        Pwet = intersect(Ps{jj}, polyshape(clipY,[zmin_all-padZ, zmin_all-padZ, z_wl, z_wl]));
        if ~isempty(Pwet.Vertices), A_sub(jj) = area(Pwet); end
    end
    qb = gamma*A_sub;
end


% derives equilibrium for a given trim 
% sum of moments and forces = 0
function R = resid_impl(v, submerged, xUnique, q_self, P, xref)
    [~,qb] = submerged(v);
    dq = qb + q_self; 
    Fres = trapz(xUnique,dq) + sum(P(:,2)); % force
    Mres = trapz(xUnique,dq.*(xUnique - xref)) + sum(P(:,2).*(P(:,1)-xref)); %moment 
    R = [Fres; Mres];
end


%% Final Arrays representing loads 

buoyantLoad = [xUnique(:), q_buoy(:)];   % [x  q_buoy]
selfWeight  = [xUnique(:), q_self(:)];   % [x  q_self]

disp('Self Weight per cross section:');
disp(selfWeight);
disp('BuoyantLoad per cross section:');
disp(buoyantLoad);