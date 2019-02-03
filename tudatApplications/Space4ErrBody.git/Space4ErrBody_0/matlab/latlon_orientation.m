function [lat_deg,lon_deg] = latlon_orientation(lat_deg,lon_deg)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



    % Convert latitude to usual range [ -90 < lat < 90 ]
    if lat_deg < 0
        lat_deg = 360 + lat_deg;
    end
    if lat_deg > 270
        lat_deg = -(360 - lat_deg);
    elseif   lat_deg > 180 && lat_deg < 270
        lat_deg = -(lat_deg - 180);
    elseif   lat_deg > 90 && lat_deg < 180
        lat_deg = 180 - lat_deg;
    end
    
        % Convert latitude to usual range [ -180 < lon < 180 ]
    if lon_deg > 180
        lon_deg = lon_deg - 360;
    elseif lon_deg < -180
        lon_deg = 360 + lon_deg;
    end






end

