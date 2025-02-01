function [H, beta, D, gamma_value] = SelfaffineFunction(x,z,K_IC_test)

    % Assuming you have your fracture profile as a vector `z` representing the height data
    % and `x` representing the horizontal position
    
    % Calculate the power spectrum of the profile
    n = length(z); % Number of data points
    z_mean = mean(z); % Remove the mean from the profile
    z_detrended = z - z_mean; % Detrend the profile
    
    % Fourier Transform to get the power spectrum
    Z = fft(z_detrended); 
    P2 = abs(Z/n).^2; % Two-sided power spectrum
    P1 = P2(1:n/2+1); % One-sided power spectrum
    P1(2:end-1) = 2*P1(2:end-1);
    
    % Frequency axis
    f = (0:(n/2))/n; 
    
    % Plotting the power spectrum
    figure(4);
    loglog(f, P1);
    xlabel('Frequency (1/Length scale)');
    ylabel('Power');
    title('Power Spectrum of the Fracture Profile');
    grid on;
    
    % Fit the power spectrum to estimate the self-affine exponent
    % Power spectrum S(f) scales as f^(-beta), where beta = 1 + 2H (H is the Hurst exponent)
    % Estimate beta by fitting a line to the log-log plot
    
    % Linear fit in log-log scale
    log_f = log(f(2:end)); % Avoid f=0 (DC component)
    log_P1 = log(P1(2:end));
    [coeffs, S] = polyfit(log_f, log_P1, 1); % Linear fit
    
    % Extract the self-affine (Hurst) exponent H
    beta = -coeffs(1); % Slope of the line
    H = (beta - 1) / 2; % Hurst exponent
    D = 2 - H; % fractal dimension, D=d+1−H, For 1D profiles d=1
    gamma_value = log(K_IC_test) / log(1 / H);
    
    % Plot the fitted line on the power spectrum
    hold on;
    fit_line = polyval(coeffs, log_f);
    plot(f(2:end), exp(fit_line), 'r--', 'LineWidth', 2);
    legend('Power Spectrum', [sprintf('Fit (H = %.2f, ', H) sprintf('B = %.2f)', beta)]);
    hold off;

    compareProfiles(z, H);

end
