function compareProfiles(real_profile, H)
    % real_profile: Input 1D real profile (vector)
    % H: Hurst exponent for generating the synthetic profile
    
    % Length of the real profile
    L = length(real_profile);
    
    % Generate the synthetic profile using the same length and Hurst exponent
    synthetic_profile = generateSyntheticProfile(L, H);
    
    % Normalize both profiles to have zero mean and unit variance for fair comparison
    real_profile_normalized = (real_profile - mean(real_profile)) / std(real_profile);
    synthetic_profile_normalized = (synthetic_profile - mean(synthetic_profile)) / std(synthetic_profile);

    % Plot the real and synthetic profiles for comparison
    figure(5);
    subplot(2, 1, 1);
    plot(real_profile_normalized, 'b', 'LineWidth', 1.5);
    hold on;
    plot(synthetic_profile_normalized, 'r--', 'LineWidth', 1.5);
    xlabel('Position');
    ylabel('Normalized Height');
    legend('Real Profile', 'Synthetic Profile');
    title(['Real vs Synthetic Self-Affine Profile (H = ', num2str(H), ')']);
    grid on;
    hold off;

    % Compare power spectra of both profiles
    subplot(2, 1, 2);
    [f_real, P_real] = computePowerSpectrum(real_profile_normalized);
    [f_synth, P_synth] = computePowerSpectrum(synthetic_profile_normalized);
    loglog(f_real, P_real, 'b', 'LineWidth', 1.5);
    hold on;
    loglog(f_synth, P_synth, 'r--', 'LineWidth', 1.5);
    xlabel('Frequency (1/Length scale)');
    ylabel('Power');
    legend('Real Profile Power Spectrum', 'Synthetic Profile Power Spectrum');
    title('Power Spectrum Comparison');
    grid on;
    hold off;

end

function z = generateSyntheticProfile(L, H)
    % L: Length of the 1D profile
    % H: Hurst exponent
    % z: Generated 1D self-affine profile
    
    % Frequency axis: We ensure that k has the same size as L
    if mod(L,2)==0
        k = [0:(L/2), -(L/2-1):-1]; % Adjust for even and odd L to ensure it matches random_phase
    elseif mod(L,2)==1
        k = [0:(L/2), -(L/2):-1];
    end
    % Power spectrum scaling (self-affine profile follows f^(-beta) or k^(-beta))
    beta = 1 + 2 * H; % Calculate beta from H
    amplitude_spectrum = k.^(-beta/2); % Amplitude spectrum proportional to k^(-beta/2)
    
    % Avoid division by zero at k=0 (set amplitude to 0 at k=0)
    amplitude_spectrum(1) = 0;
    
    % Generate random complex numbers with the defined amplitude spectrum
    random_phase = exp(1i * 2 * pi * rand(1, L)); % Random phase of size L
    Fz = amplitude_spectrum .* random_phase;      % Multiply amplitude by random phase
    
    % Apply inverse Fourier transform to generate the 1D profile
    z = real(ifft(Fz)); % Real part of the inverse FFT
    
    % Normalize the profile to have zero mean and unit variance
    z = z - mean(z);
    z = z / std(z);
end

function [f, P1] = computePowerSpectrum(profile)
    % Compute the power spectrum of the 1D profile
    L = length(profile);
    profile_detrended = profile - mean(profile); % Detrend the profile

    % Fourier Transform to get the power spectrum
    Z = fft(profile_detrended);
    P2 = abs(Z/L).^2; % Two-sided power spectrum
    P1 = P2(1:L/2+1); % One-sided power spectrum
    P1(2:end-1) = 2*P1(2:end-1); % Adjust the power for the one-sided spectrum
    
    % Frequency axis
    f = (0:(L/2))/L; 
end
