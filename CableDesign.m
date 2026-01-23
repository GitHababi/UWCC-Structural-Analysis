close all;
clc;

% Run BuoyantLoad.m and (MomentShear.m or RaceStands.m)


%% Constants
eccentricity     = 4;    % in  : above or below critical section centroid.
maxForcePerCable = 800;  % lbs : self explanatory
Fc               = 2000; % psi : max compressive strength, tensile assumed to be 10% fc
canoeHeight      = 17.4; % in  : self explanatory
lossFactor       = 0.7;  % cable strength loss factor


%% Get Cross Section Data
centroids = zeros(N);
csas      = zeros(N);
moi       = zeros(N);
keelz     = zeros(N);

for i = 1:N
    current_x = xUnique(i);
    inner_rows = (coordinates.LongitudalPosition == current_x & matches(coordinates.Side,"Inner"));
    inner_coord_table = coordinates(inner_rows,:);
    innerCoords = inner_coord_table{:,["Offset","Height"]};
    
    outer_rows = (coordinates.LongitudalPosition == current_x & matches(coordinates.Side,"Outer"));
    outer_coord_table = coordinates(outer_rows,:);
    outerCoords = outer_coord_table{:,["Offset","Height"]};
    [centroids(i), moi(i), csas(i), keelz(i)] = generateCrossSectionProperties(innerCoords,outerCoords);
end

centroids_interp = interp1(xUnique, centroids, x_common,'linear');
moi_interp       = interp1(xUnique, moi, x_common,'linear');
csas_interp      = interp1(xUnique, csas, x_common,'linear');
keelz_interp     = interp1(xUnique, keelz, x_common,'linear');
%% Calculate height of cable from bottom.
[~, MmaxPos] = min(M);
cableHeight = centroids_interp(MmaxPos) + eccentricity;


% TODO: fix this mess.
topStress = @(F,x) -M(x).*abs(canoeHeight - centroids_interp(x))./moi_interp(x) - F./csas_interp(x) - F.*(canoeHeight - cableHeight)./moi_interp(x); 
bottomStress = @(F,x) M(x).*abs(keelz_interp(x) - centroids_interp(x))./moi_interp(x) - F./csas_interp(x) + F.*abs(cableHeight - keelz_interp(x))./moi_interp(x);

% returns stresses at given Force.
% My/I - F/A - dist*F/I
