function multiscaleAnalysis(xxxq, yyyq, zzzq)

    % Perform Multiscale Analysis by downsampling at different scales
    scales = [1, 2, 4, 8];  % Different downsampling scales
    figure;
    
    for s = 1:length(scales)
        % Downsample the surface by the current scale factor
        scale = scales(s);
        xxxq_ds = imresize(xxxq, 1/scale);
        yyyq_ds = imresize(yyyq, 1/scale);
        zzzq_ds = imresize(zzzq, 1/scale);

        % Compute the power spectrum for the downsampled surface
        [f, P1] = computePowerSpectrum(zzzq_ds);

        % Fit the power spectrum to estimate the scaling exponent beta
        log_f = log(f(2:end));  % Avoid DC component
        log_P1 = log(P1(2:end));
        coeffs = polyfit(log_f, log_P1, 1);  % Linear fit
        beta = -coeffs(1);
        H = (beta - 1) / 2;  % Hurst exponent

        % Plot the log-log power spectrum
        subplot(2, ceil(length(scales)/2), s);
        loglog(f, P1, 'LineWidth', 1.5);
        xlabel('$$F$$','Interpreter','latex','FontSize',16);
        ylabel('$$PSD$$','Interpreter','latex','FontSize',16);
        title(['$$Scale: $$', num2str(scale)],'Interpreter','latex','FontSize',16);
        grid on;

        % Annotate the plot with the Beta and H values
        text(0.1, 0.9, sprintf('H = %.2f, \\beta = %.2f', H, beta), ...
            'Units', 'normalized', 'FontSize', 12);

        % Display results in the console
        disp(['Scale ', num2str(scale), ': Beta = ', num2str(beta), ', H = ', num2str(H)]);
    end
    set(gca,'FontSize',14)

end

% Function to compute the power spectrum of the surface
function [f, P1] = computePowerSpectrum(zzzq)

    % Detrend the surface profile
    zzzq_detrended = zzzq - mean(zzzq(:));

    % Perform 2D FFT and calculate the power spectrum
    Z = fft2(zzzq_detrended);
    P2 = abs(Z).^2;  % Two-sided power spectrum
    
    % One-sided power spectrum
    P1 = P2(1:size(zzzq, 1)/2+1, :);  % Symmetric portion
    
    % Radial averaging for 1D spectrum
    P1 = mean(P1, 2);  % Averaging along rows to create 1D
    
    % Frequency axis
    f = (0:(size(P1, 1)-1)) / size(P1, 1);
end
