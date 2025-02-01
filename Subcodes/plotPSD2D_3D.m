function plotPSD2D_3D(xxxq_new, yyyq_new, zzzq_new)
    % Remove NaN values by interpolation (if needed)
    zzzq_new = fillmissing(zzzq_new, 'linear');  % Fill NaN values (optional)

    % Step 1: Apply a 2D Hann Window to the surface data
    [Ny, Nx] = size(zzzq_new);  % Get the dimensions of the grid
    hann_window_x = hann(Nx);  % Hann window for x-direction
    hann_window_y = hann(Ny);  % Hann window for y-direction
    hann_window_2d = hann_window_y * hann_window_x';  % 2D Hann window
    zzzq_windowed = zzzq_new .* hann_window_2d;  % Apply the 2D window to the surface data

    % Step 2: Perform 2D Fourier Transform
    Z = fft2(zzzq_windowed);  % 2D Fourier transform of windowed surface data
    
    % Step 3: Compute the 2D Power Spectral Density (PSD)
    P2 = abs(Z / (Nx * Ny)).^2;  % Two-sided power spectrum
    P1 = fftshift(P2);  % Shift zero frequency to center for visualization

    % Step 4: Create wave vector axes for kx and ky
    kx = linspace(-Nx/2, Nx/2-1, Nx);  % Create wave vector axis for kx
    ky = linspace(-Ny/2, Ny/2-1, Ny);  % Create wave vector axis for ky
    [KX, KY] = meshgrid(kx, ky);  % Create a meshgrid for kx and ky
    k_r = sqrt(KX.^2 + KY.^2);  % Compute the radial wave vector

    % Step 5: Compute anisotropic correlation functions
    [C_pseud_x, C_pseud_y] = computeAnisoCorrelation(P1, kx, ky);
    
    % Step 6: Compute isotropic correlation function
    [C_iso, r_iso] = computeIsotropicCorrelation(P1, KX, KY);

    % Step 7: Plot anisotropic correlation functions
    plotAnisoCorrelation(kx, ky, C_pseud_x, C_pseud_y);

    % Step 8: Plot the isotropic correlation function
    figure;
    plot(r_iso, C_iso, 'LineWidth', 1.5);
    xlabel('$$Distance\;(r)$$','Interpreter','latex','FontSize',16);
    ylabel('$$C^{iso}(r)$$','Interpreter','latex','FontSize',16);
    set(gca,'FontSize',14)

    % Step 9: Plot the 3D surface of the 2D PSD
    figure;
    surf(KX, KY, log10(P1), 'EdgeColor', 'none');
    xlabel('$$k_x$$','Interpreter','latex','FontSize',16);
    ylabel('$$k_r$$','Interpreter','latex','FontSize',16);
    zlabel('$$PSD$$','Interpreter','latex','FontSize',16);
    colorbar;
    set(gca,'FontSize',14)
    view(3);  

    % Step 10: Plot C^pseudo_r versus wave vector (k_r)
    figure;
    plot(ky(:), abs(C_pseud_y(1:length(ky(:)))), 'LineWidth', 1.5);
    xlabel('Wave Vector $$k_r$$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$$C^{pseud}_r(k_r)$$', 'Interpreter', 'latex', 'FontSize', 16);
    set(gca, 'FontSize', 14);

    % Step 11: Plot C^pseudo_x versus wave vector (k_x)
    figure;
    plot(kx(:), abs(C_pseud_x), 'LineWidth', 1.5);
    xlabel('Wave Vector $$k_x$$', 'Interpreter', 'latex', 'FontSize', 16);
    ylabel('$$C^{pseud}_x(k_x)$$', 'Interpreter', 'latex', 'FontSize', 16);
    set(gca, 'FontSize', 14);

end

% Function to compute anisotropic correlation functions
function [C_pseud_x, C_pseud_y] = computeAnisoCorrelation(P1, kx, ky)
    % Get the center of the PSD
    [Ny, Nx] = size(P1);
    mid_x = round(Nx / 2);
    mid_y = round(Ny / 2);

    % Extract PSD along the k_x axis (horizontal slice)
    psd_kx = P1(mid_y, :);  % PSD along k_x axis

    % Extract PSD along the k_y axis (vertical slice)
    psd_ky = P1(:, mid_x);  % PSD along k_y axis

    % Compute inverse Fourier transform to get correlation functions
    C_pseud_x = ifft(psd_kx, 'symmetric');  % Correlation function in x-direction
    C_pseud_y = ifft(psd_ky, 'symmetric');  % Correlation function in y-direction

    % Normalize the results
    C_pseud_x = C_pseud_x / max(C_pseud_x);  % Normalize to 1
    C_pseud_y = C_pseud_y / max(C_pseud_y);  % Normalize to 1
end

% Function to compute isotropic correlation function
function [C_iso, r_iso] = computeIsotropicCorrelation(P1, KX, KY)
    % Compute the radial distance for each point in the frequency domain
    k_r = sqrt(KX.^2 + KY.^2);
    
    % Radially average the PSD to get isotropic PSD
    k_bins = linspace(0, max(k_r(:)), 100);
    P_iso = zeros(1, length(k_bins)-1);
    
    for i = 1:length(k_bins)-1
        bin_mask = (k_r >= k_bins(i)) & (k_r < k_bins(i+1));
        P_iso(i) = mean(P1(bin_mask));
    end

    % Compute the inverse Fourier transform of the isotropic PSD
    C_iso = ifft(P_iso, 'symmetric');
    
    % Create a real-space distance vector
    dk = k_bins(2) - k_bins(1);
    r_iso = (0:length(P_iso)-1) / (dk * length(P_iso));
    
    % Normalize the correlation function
    C_iso = C_iso / max(C_iso);  % Normalize to 1
end

% Function to plot anisotropic correlation functions
function plotAnisoCorrelation(kx, ky, C_pseud_x, C_pseud_y)
    % Create corresponding distance vectors (spatial domain)
    dx = 1 / max(kx);  % Spatial resolution in x-direction
    dy = 1 / max(ky);  % Spatial resolution in y-direction
    r_x = (0:length(C_pseud_x)-1) * dx;  % Real-space distances for x-direction
    r_y = (0:length(C_pseud_y)-1) * dy;  % Real-space distances for y-direction

    % Plot the anisotropic correlation function in x-direction
    figure;
    plot(r_x, C_pseud_x, 'LineWidth', 1.5);
    xlabel('$$Distance\;(r_x)$$','Interpreter','latex','FontSize',16);
    ylabel('$$C^{pseud}_x(r_x)$$','Interpreter','latex','FontSize',16);
    % title('Anisotropic Correlation Function in x-direction');
    set(gca,'FontSize',14)

    % Plot the anisotropic correlation function in y-direction
    figure;
    plot(r_y, C_pseud_y, 'LineWidth', 1.5);
    xlabel('$$Distance\;(r_r)$$','Interpreter','latex','FontSize',16);
    ylabel('$$C^{pseud}_r(r_r)$$','Interpreter','latex','FontSize',16);
    % title('Anisotropic Correlation Function in y-direction');
    set(gca,'FontSize',14)
end
