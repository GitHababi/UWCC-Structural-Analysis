function [csa, centroid] = csa(innerCoords,outerCoords)
    %CSA Calculates cross sectional area and centroid with given inner and outer
    %coordinates.
    %   Detailed explanation goes here
    
    outside = flipud(outerCoords);
    inside  = flipud(innerCoords);
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

    csa = (sum(MatrixArea_o) - sum(MatrixArea_i)) * 2;

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
    centroid = SumAz / csa;
end

