function [Ra, Rq, Rz] = surfaceRoughnessAnalysis(surfaceProfile)
    % surfaceRoughnessAnalysis calculates surface roughness parameters: Ra, Rq, and Rz.
    %
    % Inputs:
    %   surfaceProfile - A vector containing the height values of the surface profile.
    %
    % Outputs:
    %   Ra - Average Roughness
    %   Rq - Root Mean Square Roughness
    %   Rz - Peak-to-Valley Height

    % Calculate the mean height of the surface profile
    meanHeight = mean(surfaceProfile);

    % Calculate the deviations from the mean height
    deviations = surfaceProfile - meanHeight;

    % Calculate Ra (Average Roughness)
    Ra = mean(abs(deviations));

    % Calculate Rq (Root Mean Square Roughness)
    Rq = sqrt(mean(deviations.^2));

    % Calculate Rz (Peak-to-Valley Height)
    Rz = max(surfaceProfile) - min(surfaceProfile);
end
