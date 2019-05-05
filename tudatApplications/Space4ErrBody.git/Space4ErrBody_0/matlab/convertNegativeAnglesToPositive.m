function [angles] = convertNegativeAnglesToPositive( angles )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



for i = 1:numel(angles)
    % Convert angle from negative to positive
    if angles(i,1) < 0
        angles(i,1) = angles(i) + 360;
    end
    
end

end