function D = fractalDimensionAutoBox(surfaceProfile)
    % fractalDimensionAutoBox calculates the fractal dimension (D) of a
    % given surface profile using the box-counting method with automatically determined box sizes.
    %
    % Inputs:
    %   surfaceProfile - A vector containing the height values of the surface profile.
    %
    % Outputs:
    %   D - Fractal dimension

    % Ensure the surfaceProfile is a row vector
    surfaceProfile = surfaceProfile(:)';
    
    % Determine the length of the surface profile
    profileLength = length(surfaceProfile);
    
    % Define a range of box sizes automatically
    boxSizes = logspace(log10(1 / profileLength), log10(profileLength / 2), 20);
    
    % Initialize arrays to store results
    numBoxes = zeros(size(boxSizes));
    
    % Calculate the number of boxes for each box size
    for i = 1:length(boxSizes)
        boxSize = boxSizes(i);
        % Number of boxes needed to cover the profile at the current box size
        numBoxes(i) = boxCounting(surfaceProfile, boxSize);
    end
    
    % Perform linear regression on log-log plot to find the slope (fractal dimension D)
    logBoxSizes = log(1 ./ boxSizes);
    logNumBoxes = log(numBoxes);
    
    p = polyfit(logBoxSizes, logNumBoxes, 1);
    D = abs(p(1));  % Ensure the fractal dimension is positive
    
    % Plot the log-log relationship and the fitted line for visual inspection
%     figure;
%     plot(logBoxSizes, logNumBoxes, 'o');
%     hold on;
%     plot(logBoxSizes, polyval(p, logBoxSizes), '-r');
%     xlabel('log(1 / Box Size)');
%     ylabel('log(Number of Boxes)');
%     title(['Fractal Dimension (D) = ', num2str(D)]);
%     grid on;
end

function num = boxCounting(surfaceProfile, boxSize)
    % boxCounting calculates the number of boxes needed to cover the surface
    % profile at a given box size.
    %
    % Inputs:
    %   surfaceProfile - A vector containing the height values of the surface profile.
    %   boxSize - The size of the box.
    %
    % Outputs:
    %   num - Number of boxes needed to cover the profile.

    minHeight = min(surfaceProfile);
    maxHeight = max(surfaceProfile);
    range = maxHeight - minHeight;
    
    % Normalize the surface profile to a unit interval
    normalizedProfile = (surfaceProfile - minHeight) / range;
    
    % Calculate the number of boxes
    numBoxes = ceil(1 / boxSize);
    num = 0;
    
    % Traverse each box
    for i = 0:numBoxes-1
        boxMin = i * boxSize;
        boxMax = (i + 1) * boxSize;
        
        % Check if the profile intersects with the current box
        if any(normalizedProfile >= boxMin & normalizedProfile < boxMax)
            num = num + 1;
        end
    end
end
