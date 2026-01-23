function [Centroid, I, CSA, KeelHeight] = generateCrossSectionProperties(innerCoords, outerCoords)
    outside = flipud(outerCoords);
    inside  = flipud(innerCoords);
    KeelHeight = min(outside(:,2));
    t_points_o = size(outside,1);
    t_points_i = size(inside,1);
    
    MatrixArea_o = zeros(t_points_o,1);
    MatrixArea_i = zeros(t_points_i,1);
    
    for r = 1:t_points_o-1
        y1 = abs(outside(r, 1));
        y2 = abs(outside(r+1, 1));
        z1 = outside(r, 2);
        z2 = outside(r+1, 2);
        h  = z2 - z1;
        MatrixArea_o(r+1) = h * (y1 + y2) / 2;
    end
    
    for r = 1:t_points_i-1
        y1 = abs(inside(r, 1));
        y2 = abs(inside(r+1, 1));
        z1 = inside(r, 2);
        z2 = inside(r+1, 2);
        h  = z2 - z1;
        MatrixArea_i(r+1) = h * (y1 + y2) / 2;
    end
    
    CSA_half = sum(MatrixArea_o) - sum(MatrixArea_i);
    CSA      = 2 * CSA_half;
    
    Cz_o = zeros(t_points_o,1);
    TotalABar_o = zeros(t_points_o,1);
    for r = 1:t_points_o-1
        y1 = abs(outside(r, 1));
        y2 = abs(outside(r+1, 1));
        z1 = outside(r, 2);
        z2 = outside(r+1, 2);
        h  = z2 - z1;
        cz_local = h * (y1 + 2*y2) / (3*(y1 + y2));
        Cz_o(r+1) = cz_local;
        z_strip_centroid = z1 + cz_local;
        Abar = MatrixArea_o(r+1) * z_strip_centroid;
        TotalABar_o(r+1) = Abar + TotalABar_o(r);
    end
    
    Cz_i = zeros(t_points_i,1);
    TotalABar_i = zeros(t_points_i,1);
    for r = 1:t_points_i-1
        y1 = abs(inside(r, 1));
        y2 = abs(inside(r+1, 1));
        z1 = inside(r, 2);
        z2 = inside(r+1, 2);
        h  = z2 - z1;
        cz_local = h * (y1 + 2*y2) / (3*(y1 + y2));
        Cz_i(r+1) = cz_local;
        z_strip_centroid = z1 + cz_local;
        Abar = MatrixArea_i(r+1) * z_strip_centroid;
        TotalABar_i(r+1) = Abar + TotalABar_i(r);
    end
    
    SumAz = 2 * (TotalABar_o(end) - TotalABar_i(end));
    Centroid = SumAz / CSA;
    
% Moment of Inertia (around the neutral axis)
    
    Itrap_o = zeros(t_points_o,1);
    TotalI_o = zeros(t_points_o,1);
    for r = 1:t_points_o-1
        y1 = abs(outside(r, 1));
        y2 = abs(outside(r+1, 1));
        z1 = outside(r, 2);
        z2 = outside(r+1, 2);
        h  = z2 - z1;
        Itrap_o(r+1) = (h^3 * (y2^2 + 4*y2*y1 + y1^2)) / (36*(y2 + y1));
        z_strip_centroid = z1 + Cz_o(r+1);
        I_o = Itrap_o(r+1) + MatrixArea_o(r+1) * (Centroid - z_strip_centroid)^2;
        TotalI_o(r+1) = I_o + TotalI_o(r);
    end
    
    Itrap_i = zeros(t_points_i,1);
    TotalI_i = zeros(t_points_i,1);
    for r = 1:t_points_i-1
        y1 = abs(inside(r, 1));
        y2 = abs(inside(r+1, 1));
        z1 = inside(r, 2);
        z2 = inside(r+1, 2);
        h  = z2 - z1;
        Itrap_i(r+1) = (h^3 * (y2^2 + 4*y2*y1 + y1^2)) / (36*(y2 + y1));
        z_strip_centroid = z1 + Cz_i(r+1);
        I_i = Itrap_i(r+1) + MatrixArea_i(r+1) * (Centroid - z_strip_centroid)^2;
        TotalI_i(r+1) = I_i + TotalI_i(r);
    end
    
    I = 2 * (TotalI_o(end) - TotalI_i(end));
    
end

