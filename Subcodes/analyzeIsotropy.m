function analyzeIsotropy(xxxq, yyyq, zzzq)
    % Get the size of the surface grid
    [nRows, nCols] = size(zzzq);

    % Initialize arrays to store beta, H, and D values for x- and y-directions
    beta_xz = zeros(1, nRows);
    H_xz = zeros(1, nRows);
    D_xz = zeros(1, nRows);
    
    beta_yz = zeros(1, nCols);
    H_yz = zeros(1, nCols);
    D_yz = zeros(1, nCols);
    
    % Loop over each row to compute beta, H, and D for each x-z profile
    for i = 1:nRows
        profile_xz = zzzq(i, :);  % Extract the x-z profile along the row
        [beta_xz(i), H_xz(i), D_xz(i)] = estimateBetaHurstFractal(profile_xz);
    end
    
    % Loop over each column to compute beta, H, and D for each y-z profile
    for j = 1:nCols
        profile_yz = zzzq(:, j);  % Extract the y-z profile along the column
        profile_yz = profile_yz(~isnan(profile_yz));  % Remove NaN values
        if length(profile_yz) > 1
            [beta_yz(j), H_yz(j), D_yz(j)] = estimateBetaHurstFractal(profile_yz);
        else
            beta_yz(j) = NaN;
            H_yz(j) = NaN;
            D_yz(j) = NaN;
        end
    end

    % Compute mean values for x-z and y-z directions
    mean_beta_xz = nanmean(beta_xz);
    mean_H_xz = nanmean(H_xz);
    mean_D_xz = nanmean(D_xz);
    
    mean_beta_yz = nanmean(beta_yz);
    mean_H_yz = nanmean(H_yz);
    mean_D_yz = nanmean(D_yz);
    
    % H_yz(H_yz<0.2)=0.4;

    plot(xxxq(1,:)/max(xxxq(1,:)),H_yz,'LineWidth',1)
    hold on



    % Display the results for x-z and y-z directions
    fprintf('Mean values for x-z profiles:\n');
    fprintf('Mean Beta: %.3f, Mean Hurst Exponent (H): %.3f, Mean Fractal Dimension (D): %.3f\n', ...
        mean_beta_xz, mean_H_xz, mean_D_xz);

    fprintf('\nMean values for y-z profiles:\n');
    fprintf('Mean Beta: %.3f, Mean Hurst Exponent (H): %.3f, Mean Fractal Dimension (D): %.3f\n', ...
        mean_beta_yz, mean_H_yz, mean_D_yz);

    % Compare the results
    if abs(mean_beta_xz - mean_beta_yz) < 0.05 && abs(mean_H_xz - mean_H_yz) < 0.05 && abs(mean_D_xz - mean_D_yz) < 0.05
        disp('The surface is isotropic.');
    else
        disp('The surface is anisotropic.');
    end
end

% Function to estimate beta, H, and D for a given profile
function [beta, H, D] = estimateBetaHurstFractal(profile)
    % Remove the mean and detrend the profile
    profile = profile - mean(profile);
    
    % Perform FFT to get the power spectrum
    L = length(profile);
    Y = fft(profile);
    P2 = abs(Y/L).^2;    % Two-sided spectrum
    
    % Make sure indices are integers by using floor
    half_L = floor(L/2);  % This ensures integer division for half length
    P1 = P2(1:half_L+1);    % One-sided spectrum
    P1(2:end-1) = 2*P1(2:end-1);

    % Frequency axis
    f = (0:half_L)/L;  % Ensures integer size for frequency vector
    
    % Fit the power spectrum on a log-log scale to estimate beta
    log_f = log(f(2:end));  % Avoid the DC component
    log_P1 = log(P1(2:end));
    coeffs = polyfit(log_f, log_P1, 1);  % Linear fit in log-log space
    
    % Extract beta, H, and D
    beta = -coeffs(1);       % Slope gives beta
    H = (beta - 1) / 2;      % Hurst exponent
    D = 2 - H;               % Fractal dimension for 1D profile
end
