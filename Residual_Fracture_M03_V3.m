% --------------
% FractureSurf
% --------------
% Copyright © ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2024.
% All rights reserved.
%
% This code has been developed by Mohsen Talebkeikhah
% Email: mohsen.talebkeikhah@epfl.ch
%        m.talebkeikhah@gmail.com

clc
clear
close all

%% Directory containing TIFF images
image_directory = 'D:\CT-Scan\CT-Scan Data\Core-CT-scan-M03\SlicesY_2';
image_files = dir(fullfile(image_directory, '*.tif')); 

%% Rock Properties
K_IC=0.168;

%% CT-scan resolution
CT_reso = 10;  % CT-scan resolution in micron

%% Create a VideoWriter object
output_video_file = 'M03';
video_writer = VideoWriter(output_video_file);
open(video_writer);

%% Main loop
num_figures = numel(image_files);

x_surface = [];         
y_surface = []; 
z_surface1 = [];        
z_surface2 = [];

Roughness_Results = [];

% Loop through each figure
for i = 1000:1:num_figures-685 % It has been eliminated with the CT-scan data

    disp(['i = ' num2str(i) ', ' num2str(i/numel(image_files)*100) ' % of images are processed.']);
    
    % Read the TIFF image
    image = imread(fullfile(image_directory, image_files(i).name));
    
    image = rot90(image);

    % Crop Image
    image_crop=image(50:end-130,:);

    % Normalizing the image
    normalized_image = double(image_crop) / double(intmax('uint16'));

    if mod(i,10)==0
        figure(100)
        imshow(normalized_image)
        pause(0.01);
        current_frame = getframe(100);
        current_frame.cdata=current_frame.cdata(1:450,1:1720,:);
        writeVideo(video_writer, current_frame);
        pause(0.01);
    end

    % Enhancing by Applying a median filter and remove noise
    filteredImage = medfilt2(normalized_image, [3, 3]);
    enhancedImage = adapthisteq(filteredImage,'NumTiles',[8 8],'ClipLimit',0.0005);
    enhancedImage(enhancedImage(:,:)>0.3) = 1;
    enhancedImage = (enhancedImage-1)*-1;
    enhancedImage = bwareaopen(enhancedImage, 20); 
    
    % Cleaning the image by knowing that horizontal length of region should be bigger than vertical length
    % Label the connected components in the binary image
    labeledImage = bwlabel(enhancedImage);
    % Get the properties of each connected component
    regionProps = regionprops(labeledImage, 'BoundingBox');
    % Initialize an image to store the result
    cleanedImage = false(size(enhancedImage));
    % Loop through each connected component
    for k = 1 : length(regionProps)
        % Get the bounding box for this region
        thisBoundingBox = regionProps(k).BoundingBox;

        % Extract the width (x-direction) and height (y-direction) of the bounding box
        width = thisBoundingBox(3);
        height = thisBoundingBox(4);

        % Check if the width is greater than or equal to the height
        if width*0.6 >= height
            % If so, add this region to the final image
            cleanedImage(labeledImage == k) = true;
        end
    end

    % Remove the regions with horizontal length less than N number of pixel
    % Define the minimum number of pixels in the x direction
    N = 90;  

    % Label the connected components in the binary image
    labeledImage = bwlabel(cleanedImage);
    % Get the properties of each connected component
    regionProps = regionprops(labeledImage, 'BoundingBox', 'Area');
    % Initialize an image to store the result
    finalImage = false(size(cleanedImage));
    % Loop through each connected component
    for k = 1 : length(regionProps)
        % Get the bounding box for this region
        thisBoundingBox = regionProps(k).BoundingBox;

        % Check if the width (x direction) of the bounding box is greater than N
        if thisBoundingBox(3) >= N
            % If so, add this region to the final image
            finalImage(labeledImage == k) = true;
        end
    end

    % Apply Canny edge detection
    edges = edge(finalImage, 'Canny');
    
    % Find the row and column indices of the non-zero elements (edges)
    [row_indices, col_indices] = find(edges);
    
    if size(row_indices,1)<1
        continue;
    end
    
    % Plot the processed image within the loop
    figure(1)
    subplot(3,2,1)
    imshow(normalized_image) 
    title('normalized image')
    subplot(3,2,3)
    imshow(enhancedImage)
    title('enhanced image')
    subplot(3,2,5)
    imshow(cleanedImage)
    title('cleaned image')
    subplot(3,2,2)
    imshow(finalImage)
    title('final image')
    subplot(3,2,4)
    imshow(edges);
    title('detected edges');
    subplot(3,2,6)
    imshow(normalized_image);
    title('normalized image + detected edges');
    hold on;
    plot(col_indices, row_indices, 'r.','markersize',1); 
    hold off;
    pause(0.1)
    
    lin1x=[];
    lin1y=[];
    lin2x=[];
    lin2y=[];
    midlinx=[];
    midliny=[];
    zz=1;
    for kk=1:size(finalImage,2)
        ind=find(col_indices==kk);
        if numel(ind)>1 
            if abs(row_indices(max(ind))-row_indices(min(ind)))<(2*CT_reso)
                lin1x(zz)=kk*CT_reso;
                lin1y(zz)=row_indices(min(ind))*CT_reso;
                lin2x(zz)=kk*CT_reso;
                lin2y(zz)=row_indices(max(ind))*CT_reso;
                midlinx(zz)=kk*CT_reso;
                midliny(zz)=(lin1y(zz)+lin2y(zz))/2;
                zz=zz+1;
            end
        end
    end

    % Store the x, y, and z values for this loop - mm
    x_surface = [x_surface, lin1x/1000];
    y_surface = [y_surface, i * CT_reso * ones(1, length(lin1x))/1000]; 

    % 1st surface - mm
    z_surface1 = [z_surface1, lin1y/1000];
    
    % 2nd surface - mm       
    z_surface2 = [z_surface2, lin2y/1000];

    % Calculate the autocorrelation of z
    [acf, lags] = autocorr(midliny,floor(numel(midliny)*0.8));
    
    % Plot the autocorrelation
    figure(2)
    plot(lags, acf, 'k', 'LineWidth', 1.5);
    hold on;

    p_acf_lag = polyfit(lags(1:floor(numel(midliny)*0.05)),acf(1:floor(numel(midliny)*0.05)),1);
    x_acf = lags(1):lags(floor(numel(midliny)*0.2));
    y_acf = polyval(p_acf_lag,x_acf);

    plot(x_acf, y_acf, 'r--', 'LineWidth', 1.5);
    title(['slope = ' num2str(p_acf_lag(1))]);
    xlabel('Lag');
    ylabel('Autocorrelation');
    grid on;
    hold off;

    Pearson_Corr = corr(lin1y', lin2y');

    % Opening perpendicular to the mid-surface
    W_pr = NormDistFunc([midlinx',midliny'],[lin1x',lin1y'])+...
           NormDistFunc([midlinx',midliny'],[lin2x',lin2y']);

    W_pr_nonan=W_pr(~isnan(W_pr));
    
    figure(3)
    subplot(3,1,1)
    plot(lin1x,lin1y,'r.','markersize',3)
    hold on
    plot(lin2x,lin2y,'b.','markersize',3)
    plot(midlinx,midliny,'k.','markersize',3)
    xlim([0 16020])
    title('detected fracture surface')
    hold off;

    subplot(3,1,2)
    plot(lin1x,W_pr,'k.','markersize',3);
    yline(median(W_pr_nonan),'m--')
    xlim([0 16020])
    title('fracture opening');
    hold off;

    % Calculate mean and standard deviation
    data_mean = mean(W_pr_nonan);
    data_std = std(W_pr_nonan);
    data_med = median(W_pr_nonan);
    
    subplot(3,1,3)
    [hist_counts, hist_edges] = histcounts(W_pr_nonan, 30); 
    histogram(W_pr_nonan, 30, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'k'); 
    bin_centers = hist_edges(1:end-1) + diff(hist_edges)/2; % Calculate bin centers
    [max_count, max_idx] = max(hist_counts); % Find the peak count and its index
    crest_value = bin_centers(max_idx); % Get the corresponding bin center

    hold on;
    
    % Add mean and standard deviation lines
    y_limits = ylim; % Get current y-axis limits
    line([data_mean, data_mean], y_limits, 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Mean');
    line([crest_value, crest_value], y_limits, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Crest');
    line([data_med, data_med], y_limits, 'Color', 'm', 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Median');
    line([data_mean - data_std, data_mean - data_std], y_limits, 'Color', 'g', 'LineWidth', 1.5, 'LineStyle', ':', 'DisplayName', '-1 Std Dev');
    line([data_mean + data_std, data_mean + data_std], y_limits, 'Color', 'g', 'LineWidth', 1.5, 'LineStyle', ':', 'DisplayName', '+1 Std Dev');
    
    title(sprintf('Histogram with Mean = %.2f and Std Dev = %.2f', data_mean, data_std));
    xlabel('Data Values');
    ylabel('Frequency');
    grid on;
    legend('Histogram', 'Mean', 'Crest', '-1 Std Dev', '+1 Std Dev', 'Location', 'best');
    hold off;
    pause(0.1)
    
    % Roughness Analysis
    % Self-Affinement - Power Spectrum Method
    [H1, beta1, D1, gama1] = SelfaffineFunction(lin1x,lin1y,K_IC);

    [H2, beta2, D2, gama2] = SelfaffineFunction(lin2x,lin2y,K_IC);

    % Roughness properties
    % mid-surf
    [Ra_v, Rq_v, Rz_v] = SurfaceRoughnessAnalysis(midliny);
    
    % surf 1
    D_m1 = NormDistFunc([midlinx',midliny'],[lin1x',lin1y']);
    D_m1 = D_m1(~isnan(D_m1));

    [Ra_m1, Rq_m1, Rz_m1] = SurfaceRoughnessAnalysis(D_m1);
    
    % surf 2
    D_m2 = NormDistFunc([midlinx',midliny'],[lin2x',lin2y']);
    D_m2 = D_m2(~isnan(D_m2));

    [Ra_m2, Rq_m2, Rz_m2] = SurfaceRoughnessAnalysis(D_m2);

    Roughness_Results = [Roughness_Results; [i*CT_reso, data_med, crest_value, data_mean, data_std, p_acf_lag(1), Pearson_Corr, ...
                                             H1, beta1, D1, real(gama1), H2, beta2, D2, real(gama1),...
                                             Ra_v, Rq_v, Rz_v, Ra_m1, Rq_m1, Rz_m1, Ra_m2, Rq_m2, Rz_m2]];
end

% Close the video file
close(video_writer);

%% Export the scattered data & Save the roughness data
% save the real position of cloud points
writematrix([x_surface', ((y_surface-y_surface(end))*(-1))', z_surface1'],'M03_Fracture_Surfaces_Upper_real_position.csv');
writematrix([x_surface', ((y_surface-y_surface(end))*(-1))', z_surface2'],'M03_Fracture_Surfaces_Lower_real_position.csv');
writematrix(Roughness_Results,'M03_Roughness_Results.csv');

%% Load data from .csv
data_lower = csvread('M03_Fracture_Surfaces_Lower_real_position.csv');
data_upper = csvread('M03_Fracture_Surfaces_Upper_real_position.csv');
Roughness_Data = csvread('M03_Roughness_Results.csv');

x_surface=data_upper(:,1);
y_surface1=data_upper(:,2);
z_surface1=data_upper(:,3);
y_surface2=data_lower(:,2);
z_surface2=data_lower(:,3);

%% H, std, corrolation coefficients, and residual opening
clc
close all

X_Corr = ((Roughness_Data(:,1)-Roughness_Data(end,1))*-1)/(10^(4));
X_Corr_norm = X_Corr/max(X_Corr);

figure
plot(X_Corr_norm,Roughness_Data(:,8),'b','DisplayName','upper profile');
hold on;
plot(X_Corr_norm,Roughness_Data(:,12),'r','DisplayName','lower profile');
% xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
% xlabel('$$r/r_{max}$$','Interpreter','latex','FontSize',16);
ylabel('$$H$$','Interpreter','latex','FontSize',16)
legend show;
set(gca,'FontSize',14)

figure
plot(X_Corr,Roughness_Data(:,5),'b');
hold on;
plot(X_Corr,smooth(Roughness_Data(:,5),0.1),'k');
xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
ylabel('$$\sigma$$','Interpreter','latex','FontSize',16)
set(gca,'FontSize',14)

% Corrolation Coefficients
figure
plot(X_Corr,Roughness_Data(:,6),'b');
hold on;
plot(X_Corr,smooth(Roughness_Data(:,6),0.01),'k', 'LineWidth', 1.5);
xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
ylabel('$$r_{acf}$$','Interpreter','latex','FontSize',16)
set(gca,'FontSize',14)

figure
plot(X_Corr,Roughness_Data(:,7),'b');
hold on;
plot(X_Corr,smooth(Roughness_Data(:,7),0.01),'k', 'LineWidth', 1.5);
xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
ylabel('$$R^{2}$$','Interpreter','latex','FontSize',16)
set(gca,'FontSize',14)

figure
yyaxis left   
hold on;
plot(X_Corr_norm,smooth(Roughness_Data(:,7),0.1), 'LineWidth', 1.5);
% xlabel('$$r/r_{max}$$','Interpreter','latex','FontSize',16);
ylabel('$$R_{P}$$','Interpreter','latex','FontSize',16)
set(gca,'FontSize',16)

yyaxis right   
hold on;
plot(X_Corr_norm,smooth(Roughness_Data(:,6),0.1), 'LineWidth', 1.5);
% xlabel('$$r/r_{max}$$','Interpreter','latex','FontSize',16);
ylabel('$$R_{acf}$$','Interpreter','latex','FontSize',16)
set(gca,'FontSize',16)

% Residual opening
figure
plot(X_Corr,Roughness_Data(:,2),'b');
hold on;
plot(X_Corr,smooth(Roughness_Data(:,2),0.1),'k', 'LineWidth', 1.5);
xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
ylabel('$$w\;(\mu m)$$','Interpreter','latex','FontSize',16)
set(gca,'FontSize',16)

Data_Pro=[X_Corr, smooth(X_Corr,Roughness_Data(:,2),0.5)];

error=CT_reso*ones(1,size(Data_Pro,1));
figure
e=errorbar(Data_Pro(:,1), Data_Pro(:,2), error, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
e.Color = [220 220 220]/255;
hold on
plot(Data_Pro(:,1), Data_Pro(:,2),'linewidth',1.5,'Color','k');
xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
ylabel('$$w\;(\mu m)$$','Interpreter','latex','FontSize',16)
set(gca,'FontSize',16)

% General roughness analysis
figure
plot(X_Corr,Roughness_Data(:,16))
hold on
plot(X_Corr,Roughness_Data(:,17))
plot(X_Corr,Roughness_Data(:,18))
legend('Average Roughness','Root Mean Square Roughness','Peak-to-Valley Height')
ylabel('$$\alpha\;(\mu m)$$','Interpreter','latex','FontSize',16)
xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
title('mid-surface')
set(gca,'FontSize',14)

figure
plot(X_Corr,Roughness_Data(:,19))
hold on
plot(X_Corr,Roughness_Data(:,20))
plot(X_Corr,Roughness_Data(:,21))
legend('Average Roughness','Root Mean Square Roughness','Peak-to-Valley Height')
xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
ylabel('$$\alpha\;(\mu m)$$','Interpreter','latex','FontSize',16)
title('upper surface')
set(gca,'FontSize',14)

figure
plot(X_Corr,Roughness_Data(:,22))
hold on
plot(X_Corr,Roughness_Data(:,23))
plot(X_Corr,Roughness_Data(:,24))
legend('Average Roughness','Root Mean Square Roughness','Peak-to-Valley Height')
xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
ylabel('$$\alpha\;(\mu m)$$','Interpreter','latex','FontSize',16)
title('lower surface')
set(gca,'FontSize',14)


%% 3D fracture suface
clc
close all

% 1st profile - cloud points
% figure
% scatter3(x_surface, y_surface1, z_surface1, 10, z_surface1, 'filled');
% xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
% ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
% zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
% title('upper profile');
% % colorbar;
% ax = gca; 
% ax.DataAspectRatio = [0.8 0.5 0.3]; 

% 1st profile - surface
% Your original data
xxx1 = x_surface;
yyy1 = y_surface1;
zzz1 = z_surface1;

% Create a meshgrid
[xxxq1, yyyq1] = meshgrid(min(xxx1):0.1:max(xxx1), min(yyy1):0.1:max(yyy1));
zzzq1 = griddata(xxx1, yyy1, zzz1, xxxq1, yyyq1, 'cubic');

% Mask NaN values for the surface plot
% figure
% surf(xxxq1, yyyq1, zzzq1, 'EdgeColor', 'none'); % No edges, smoother plot
% 
% % Set NaN points as transparent
% set(gca, 'ALim', [min(zzzq1(:)), max(zzzq1(:))]);  % Set color limits
% alpha(surf(xxxq1, yyyq1, zzzq1), double(~isnan(zzzq1))); % Make NaN points transparent
% 
% % Set labels and title
% xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
% ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
% zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
% title('upper profile');
% 
% % Adjust axis properties
% ax = gca; 
% ax.DataAspectRatio = [0.8 0.5 0.3]; 

% 1st profile - without grid lines
% Plot surface without gridlines
figure
surf(xxxq1, yyyq1, zzzq1, 'EdgeColor', 'none');  % 'EdgeColor', 'none' removes the gridlines

% Set NaN points as transparent
alpha(surf(xxxq1, yyyq1, zzzq1, 'EdgeColor', 'none'), double(~isnan(zzzq1)));  % Transparent NaN values

% Set labels and title
xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
title('upper profile');

% Adjust axis properties
ax = gca; 
ax.DataAspectRatio = [0.8 0.5 0.3]; 


% 2st profile - cloud points
% figure
% scatter3(x_surface, y_surface2, z_surface2, 10, z_surface2, 'filled');
% xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
% ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
% zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
% title('lower profile');
% % colorbar;
% ax = gca; 
% ax.DataAspectRatio = [0.8 0.5 0.3]; 

% 2nd profile - surface
% Your original data
xxx2 = x_surface;
yyy2 = y_surface2;
zzz2 = z_surface2;

% Create a meshgrid
[xxxq2, yyyq2] = meshgrid(min(xxx2):0.1:max(xxx2), min(yyy2):0.1:max(yyy2));
zzzq2 = griddata(xxx2, yyy2, zzz2, xxxq2, yyyq2, 'cubic');

% Mask NaN values for the surface plot
% figure
% surf(xxxq2, yyyq2, zzzq2, 'EdgeColor', 'none'); % No edges, smoother plot
% 
% % Set NaN points as transparent
% set(gca, 'ALim', [min(zzzq2(:)), max(zzzq2(:))]);  % Set color limits
% alpha(surf(xxxq2, yyyq2, zzzq2), double(~isnan(zzzq2))); % Make NaN points transparent
% 
% % Set labels and title
% xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
% ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
% zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
% title('lower profile');
% 
% % Adjust axis properties
% ax = gca; 
% ax.DataAspectRatio = [0.8 0.5 0.3]; 

% 2nd profile - without grid lines
% Plot surface without gridlines
figure
surf(xxxq2, yyyq2, zzzq2, 'EdgeColor', 'none');  % 'EdgeColor', 'none' removes the gridlines

% Set NaN points as transparent
alpha(surf(xxxq2, yyyq2, zzzq2, 'EdgeColor', 'none'), double(~isnan(zzzq2)));  % Transparent NaN values

% Set labels and title
xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
title('lower profile');

% Adjust axis properties
ax = gca; 
ax.DataAspectRatio = [0.8 0.5 0.3]; 

%% two surfaces together
clc
close all;

offset = 2; % mm

figure
colormap(parula);
h1 = surf(xxxq1, yyyq1, zzzq1 + offset, 'EdgeColor', 'none');  
alpha(h1, double(~isnan(zzzq1))); 
alpha(h1,0.7)

hold on;

colormap(jet);
h2 = surf(xxxq2, yyyq2, zzzq2, 'EdgeColor', 'none');  
alpha(h2, double(~isnan(zzzq2)));  


xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);

% title('lower profile');

ax = gca; 
ax.DataAspectRatio = [0.8 0.5 0.3]; 
% colorbar

%% Contact points in CT-scan
clc
close all

epsilon = 0.02;

% Calculate the height difference between the two profiles
delta_z = zzz2 - zzz1;

y_connect = yyy1(delta_z<=epsilon);
x_connect = xxx1(delta_z<=epsilon);

y_connect_coord = Roughness_Data(end,1)/CT_reso - y_connect*1000/CT_reso;
x_connect_coord = x_connect*1000/CT_reso;

num_con_pon = 6000; % 6000, 4000, 2000, 1000, 999, 400, 1
index_CT_image = floor(y_connect_coord(num_con_pon));

image = imread(fullfile(image_directory, image_files(index_CT_image).name));
image = rot90(image);
image_crop = image(50:end-130,:);

figure
imshow(image_crop)
hold on; % Keep the image for adding more plots
% line([x_connect_coord(num_con_pon) x_connect_coord(num_con_pon)], [1 size(image_crop,1)], 'Color', 'r'); % Draw the vertical line

% Plot the scale
[imgHeight, imgWidth, numChannels] = size(image_crop);

scaleBarLength = 100;  % 1 mm = 1000 microns -> 100 pixels
scaleBarHeight = 5;    % Height of the bar in pixels
xPosition = 20;        % X-coordinate for the scale bar (left padding)
yPosition = 60;        % Y-coordinate (bottom padding)

% Draw the white rectangle
imgWithBar = insertShape(image_crop, 'FilledRectangle', [xPosition, yPosition, scaleBarLength, scaleBarHeight], ...
                         'Color', 'white', 'Opacity', 1);

imshow(imgWithBar);


%% Calculation of contact area - full fracture surface
clc
close all

delta_z = zzzq2 - zzzq1;

for i123 = 1:size(delta_z,1)
    for j123 = 1:size(delta_z,2)
        if delta_z(i123,j123)>0.3 | delta_z(i123,j123)<0
            delta_z(i123,j123) = nan;
        end
    end
end

% Plot the contact points on the surface profiles' difference
figure;
colormap(parula);
h1 = surf(xxxq1, yyyq1, delta_z, 'EdgeColor', 'none'); 
alpha(h1, double(~isnan(delta_z)));  % Transparent NaN values
hold on;

xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
ax = gca; 
ax.DataAspectRatio = [0.6 0.8 0.05]; 
view(2)


epsilon = 0.02;  
contact_points = delta_z <= epsilon;
dx = xxxq1(1, 2) - xxxq1(1, 1);  % Grid spacing in x-direction
dy = yyyq1(2, 1) - yyyq1(1, 1);  % Grid spacing in y-direction
area_per_point = dx * dy;        % Area of each grid cell


A_0 = sum(contact_points(:)) * area_per_point;  % Total contact area
A_total = sum(sum(ones(size(delta_z)))) * area_per_point; % Total fracture area
fprintf('Real Contact Area: %.3f mm^2\n', A_0);
fprintf('Real Total Area: %.3f mm^2\n', A_total);


scatter3(xxxq1(contact_points), yyyq1(contact_points), zzzq1(contact_points), ...
    3, 'filled', 'MarkerFaceColor', 'r');
colorbar;

% Plot the contact points on the surface profiles
Pixel_threshold = 617;
for i12=1:size(delta_z,1)
    for j12=1:size(delta_z,2)
            xxxq11(i12,j12)=xxxq1(i12,j12);
            yyyq11(i12,j12)=yyyq1(i12,j12);
            if i12>Pixel_threshold & isnan(delta_z(i12,j12))
                delta_z(i12,j12)=0.019;
                zzzq11(i12,j12)=2.8;
            else
                zzzq11(i12,j12)=zzzq1(i12,j12);
            end
    end
end

figure;
surf(xxxq1, yyyq1, zzzq1, 'EdgeColor', 'none'); 
hold on;
surf(xxxq2, yyyq2, zzzq2, 'EdgeColor', 'none');
xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
% colorbar;
ax = gca; 
ax.DataAspectRatio = [0.8 0.5 1]; 
scatter3(xxxq11(contact_points), yyyq11(contact_points), zzzq11(contact_points), ...
    3, 'filled', 'MarkerFaceColor', 'r');


% Calculate the number of contact points along the y-axis
num_contact_points_y = sum(contact_points, 2);  % Sum across each row (y-direction)

% Plot the number of contact points versus the y-axis
figure;
plot(yyyq1(:,1)/max(yyyq1(:,1)), num_contact_points_y/size(yyyq1,2), 'LineWidth', 1.5);
xlabel('$$r/r_{max}$$','Interpreter','latex','FontSize',16);
ylabel('$$P_{c}/P_{tx}$$','Interpreter','latex','FontSize',16);
set(gca,'FontSize',14)

figure;
plot(yyyq1(:,1)/max(yyyq1(:,1)), num_contact_points_y/size(yyyq1,2), 'LineWidth', 1.5,'color',[0.8, 0.8, 0.8]);
hold on;
plot(yyyq1(:,1)/max(yyyq1(:,1)), movmean(num_contact_points_y/size(yyyq1,2),5), 'LineWidth', 1.5,'color','k');
xlabel('$$r/r_{max}$$','Interpreter','latex','FontSize',16);
ylabel('$$P_{c}/P_{tx}$$','Interpreter','latex','FontSize',16);
set(gca,'FontSize',14)

figure;
plot(yyyq1(:,1), num_contact_points_y, 'LineWidth', 1.5);
xlabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
ylabel('$$Number\;of\;Contact\;Points$$','Interpreter','latex','FontSize',16);
set(gca,'FontSize',14)

%% Calculation of contact area - Covering the nan parts at the fracture tip
clc
close all

delta_z = zzzq2 - zzzq1;

Pixel_threshold = 617;
for i12=1:size(delta_z,1)
    for j12=1:size(delta_z,2)
            xxxq11(i12,j12)=xxxq1(i12,j12);
            yyyq11(i12,j12)=yyyq1(i12,j12);
            if i12>Pixel_threshold & isnan(delta_z(i12,j12))
                delta_z(i12,j12)=0.019;
                zzzq11(i12,j12)=2.8;
            else
                zzzq11(i12,j12)=zzzq1(i12,j12);
            end
    end
end

epsilon = 0.02;  
contact_points = delta_z <= epsilon;
dx = xxxq1(1, 2) - xxxq1(1, 1);  % Grid spacing in x-direction
dy = yyyq1(2, 1) - yyyq1(1, 1);  % Grid spacing in y-direction
area_per_point = dx * dy;  % Area of each grid cell
A_0 = sum(contact_points(:)) * area_per_point;  % Total contact area
A_total = sum(sum(ones(size(delta_z)))) * area_per_point; % Total fracture area
fprintf('Real Contact Area: %.3f mm^2\n', A_0);
fprintf('Real Total Area: %.3f mm^2\n', A_total);

% Plot the contact points on the surface profiles
figure;
surf(xxxq1, yyyq1, zzzq1, 'EdgeColor', 'none'); 
hold on;
surf(xxxq2, yyyq2, zzzq2, 'EdgeColor', 'none');
xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
% colorbar;
ax = gca; 
ax.DataAspectRatio = [0.8 0.5 1]; 

scatter3(xxxq11(contact_points), yyyq11(contact_points), zzzq11(contact_points), ...
    3, 'filled', 'MarkerFaceColor', 'r');


figure;
surf(xxxq1, yyyq1, zzzq1, 'EdgeColor', 'none'); 
hold on;
surf(xxxq2, yyyq2, zzzq2, 'EdgeColor', 'none');
xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
% colorbar;
ax = gca; 
ax.DataAspectRatio = [0.8 0.5 1]; 
scatter3(xxxq1(contact_points), yyyq1(contact_points), zzzq1(contact_points), ...
    3, 'filled', 'MarkerFaceColor', 'r');

num_contact_points_y = sum(contact_points, 2);  % Sum across each row (y-direction)

figure;
plot(yyyq1(:,1)/max(yyyq1(:,1)), num_contact_points_y/size(yyyq1,2), 'LineWidth', 1.5);
xlabel('$$r/r_{max}$$','Interpreter','latex','FontSize',16);
ylabel('$$P_{c}/P_{tx}$$','Interpreter','latex','FontSize',16);
set(gca,'FontSize',14)

figure;
plot(yyyq1(:,1), num_contact_points_y, 'LineWidth', 1.5);
xlabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
ylabel('$$Number\;of\;Contact\;Points$$','Interpreter','latex','FontSize',16);
set(gca,'FontSize',14)

%% Define a specific rectangular area which cover the most of the surface 
clc
close all

% profile 1
figure
surf(xxxq1, yyyq1, zzzq1, 'EdgeColor', 'none');  % 'EdgeColor', 'none' removes the gridlines

% Set NaN points as transparent
alpha(surf(xxxq1, yyyq1, zzzq1, 'EdgeColor', 'none'), double(~isnan(zzzq1)));  % Transparent NaN values

% Set labels and title
xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
title('upper profile');

% Adjust axis properties
ax = gca; 
ax.DataAspectRatio = [0.8 0.5 0.3]; 

% profile 2
figure
surf(xxxq2, yyyq2, zzzq2, 'EdgeColor', 'none');  % 'EdgeColor', 'none' removes the gridlines

% Set NaN points as transparent
alpha(surf(xxxq2, yyyq2, zzzq2, 'EdgeColor', 'none'), double(~isnan(zzzq2)));  % Transparent NaN values

% Set labels and title
xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
title('lower profile');

% Adjust axis properties
ax = gca; 
ax.DataAspectRatio = [0.8 0.5 0.3]; 

% Define the region of interest
x_min = 0.32;
x_max = 15.22;
y_min = 0.5;
y_max = 62.1;

% Find indices where xxxq and yyyq are within the specified range
valid_x_idx1 = xxxq1(1,:) >= x_min & xxxq1(1,:) <= x_max;  % Indices for x-direction
valid_y_idx1 = yyyq1(:,1) >= y_min & yyyq1(:,1) <= y_max;  % Indices for y-direction

valid_x_idx2 = xxxq2(1,:) >= x_min & xxxq2(1,:) <= x_max;  % Indices for x-direction
valid_y_idx2 = yyyq2(:,1) >= y_min & yyyq2(:,1) <= y_max;  % Indices for y-direction

% Extract the corresponding submatrices for xxxq, yyyq, and zzzq
xxxq_new1 = xxxq1(valid_y_idx1, valid_x_idx1);
yyyq_new1 = yyyq1(valid_y_idx1, valid_x_idx1);
zzzq_new1 = zzzq1(valid_y_idx1, valid_x_idx1);

xxxq_new2 = xxxq2(valid_y_idx2, valid_x_idx2);
yyyq_new2 = yyyq2(valid_y_idx2, valid_x_idx2);
zzzq_new2 = zzzq2(valid_y_idx2, valid_x_idx2);

% Plot the result
figure
surf(xxxq_new1, yyyq_new1, zzzq_new1, 'EdgeColor', 'none');  
xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
title('Smaller Meshgrid without NaNs - profile 1');
ax = gca; 
ax.DataAspectRatio = [0.8 0.5 1]; 

figure
surf(xxxq_new2, yyyq_new2, zzzq_new2, 'EdgeColor', 'none'); 
xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
title('Smaller Meshgrid without NaNs - profile 2');
ax = gca; 
ax.DataAspectRatio = [0.8 0.5 1]; 

%% Calculation of contact area - rectangular area which cover the most of the surface
% clc
% close all
% 
% delta_z = zzzq_new2 - zzzq_new1;
% 
% epsilon = 0.02;
% contact_points = delta_z <= epsilon;
% dx = xxxq_new1(1, 2) - xxxq_new1(1, 1);  % Grid spacing in x-direction
% dy = yyyq_new1(2, 1) - yyyq_new1(1, 1);  % Grid spacing in y-direction
% area_per_point = dx * dy;  % Area of each grid cell
% A_0 = sum(contact_points(:)) * area_per_point;  % Total contact area
% 
% % Display the real contact area
% fprintf('Real Contact Area: %.3f mm^2\n', A_0);
% 
% figure;
% surf(xxxq_new1, yyyq_new1, zzzq_new1, 'EdgeColor', 'none'); 
% hold on;
% surf(xxxq_new2, yyyq_new2, zzzq_new2, 'EdgeColor', 'none');
% xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
% ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
% zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
% % colorbar;
% ax = gca; 
% ax.DataAspectRatio = [0.8 0.5 1]; 
% 
% scatter3(xxxq_new1(contact_points), yyyq_new1(contact_points), zzzq_new1(contact_points), ...
%     3, 'filled', 'MarkerFaceColor', 'r');
% 
% % Calculate the number of contact points along the y-axis
% num_contact_points_y = sum(contact_points, 2);  % Sum across each row (y-direction)
% 
% % Plot the number of contact points versus the y-axis
% figure;
% plot(yyyq_new1(:,1)/max(yyyq_new1(:,1)), num_contact_points_y/size(yyyq_new1,2), 'LineWidth', 1.5);
% xlabel('$$r/r_{max}\;(mm)$$','Interpreter','latex','FontSize',16);
% ylabel('$$P_{c}/P_{tx}$$','Interpreter','latex','FontSize',16);
% set(gca,'FontSize',14)
% 
% figure;
% plot(yyyq_new1(:,1), num_contact_points_y, 'LineWidth', 1.5);
% xlabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
% ylabel('$$Number\;of\;Contact\;Points$$','Interpreter','latex','FontSize',16);
% set(gca,'FontSize',14)
% 
% figure;
% plot(yyyq_new1(:,1)/max(yyyq_new1(:,1)), num_contact_points_y/size(yyyq_new1,2), 'LineWidth', 1.5,'color',[0.8, 0.8, 0.8]);
% hold on;
% plot(yyyq_new1(:,1)/max(yyyq_new1(:,1)), movmean(num_contact_points_y/size(yyyq_new1,2),5), 'LineWidth', 1.5,'color','k');
% xlabel('$$r/r_{max}\;(mm)$$','Interpreter','latex','FontSize',16);
% ylabel('$$P_{c}/P_{tx}$$','Interpreter','latex','FontSize',16);
% set(gca,'FontSize',14)

%% h_rms, h'_rms, and h''_rms with surface data
clc
close all

% Assuming you have the grid data for the surface profile:
% xxxq_new1, yyyq_new1: X and Y coordinates of the grid points
% zzzq_new1: Z values (height) of the surface profile

% Calculate h_rms (root-mean-square height)
mean_z = mean(zzzq_new2(:));  % Mean height
h_rms = sqrt(mean((zzzq_new2(:) - mean_z).^2));  % RMS height

% Calculate h'_rms (root-mean-square slope)
% First, compute the derivatives in the x and y directions
[dz_dx, dz_dy] = gradient(zzzq_new2, mean(diff(xxxq_new2(1,:))), mean(diff(yyyq_new2(:,1))));

% Calculate the RMS slope
h_rms_prime = sqrt(mean(dz_dx(:).^2 + dz_dy(:).^2));  % RMS slope

% Calculate h''_rms (root-mean-square curvature)
% Second derivatives of the surface in the x and y directions
[d2z_dx2, d2z_dxdy] = gradient(dz_dx, mean(diff(xxxq_new2(1,:))), mean(diff(yyyq_new2(:,1))));
[~, d2z_dy2] = gradient(dz_dy, mean(diff(xxxq_new2(1,:))), mean(diff(yyyq_new2(:,1))));

% Calculate the RMS curvature
h_rms_double_prime = sqrt(mean(d2z_dx2(:).^2 + d2z_dy2(:).^2 + 2 * d2z_dxdy(:).^2));  % RMS curvature

% Display the results
fprintf('h_rms (RMS height): %.4f\n', h_rms);
fprintf('h''_rms (RMS slope): %.4f\n', h_rms_prime);
fprintf('h''''_rms (RMS curvature): %.4f\n', h_rms_double_prime);

%% h_rms, h'_rms, and h''_rms with PSD

close all;

% Assuming you have the grid data for the surface profile:
% xxxq_new1, yyyq_new1: X and Y coordinates of the grid points in mm
% zzzq_new1: Z values (height) of the surface profile in mm

% Step 1: Define grid spacing in the X and Y directions
dx = mean(diff(xxxq_new1(1, :)));  % Grid spacing in x-direction in mm
dy = mean(diff(yyyq_new1(:, 1)));  % Grid spacing in y-direction in mm

% Step 2: Compute the 2D Fourier transform of the surface heights
fft2_z = fft2(zzzq_new1);  % 2D Fourier transform of the height data

[nRows, nCols] = size(zzzq_new1);

% Compute the Power Spectral Density (PSD)
psd2_z = (abs(fft2_z).^2) / (nRows * nCols);  % Normalize by the number of points

% Multiply by grid spacing (dx, dy) to account for physical dimensions (mm)
psd2_z = psd2_z * (dx * dy);  % Now the PSD has units of mm^2

% Normalize by the total area of the grid (to adjust units)
A = (nRows * dx) * (nCols * dy);  % Total area in mm^2
psd2_z = psd2_z / A;  % PSD now has the correct units: mm^2/mm^2 = dimensionless

% Handle odd/even lengths for frequency grids
if mod(nCols, 2) == 0
    kx = (2 * pi / (nCols * dx)) * [0:(nCols/2 - 1), -nCols/2:-1];
else
    kx = (2 * pi / (nCols * dx)) * [0:(nCols-1)/2, -(nCols-1)/2:-1];
end

if mod(nRows, 2) == 0
    ky = (2 * pi / (nRows * dy)) * [0:(nRows/2 - 1), -nRows/2:-1];
else
    ky = (2 * pi / (nRows * dy)) * [0:(nRows-1)/2, -(nRows-1)/2:-1];
end

% Create a 2D meshgrid for frequencies
[KX, KY] = meshgrid(kx, ky);

% Compute the radial frequency (wavevector magnitude q) in mm^-1
q = sqrt(KX.^2 + KY.^2);
q = fftshift(q);  % Align with the PSD

% Step 4: Calculate the integrals for h_rms, h'_rms, and h''_rms
% Calculate h_rms (RMS height in mm)
h_rms = sqrt(sum(psd2_z(:)));  % Summing over all frequencies, result in mm

% Calculate h'_rms (RMS slope dimensionless)
h_rms_prime = sqrt(sum(psd2_z(:)));  % Dimensionless because no q is involved

% Calculate h''_rms (RMS curvature in mm^-1)
h_rms_double_prime = sqrt(sum((q(:).^2) .* psd2_z(:)));  % Units of mm^-1

% Display the results
fprintf('h_rms (RMS height in mm): %.4f mm\n', h_rms);
fprintf('h''_rms (RMS slope dimensionless): %.4f\n', h_rms_prime);
fprintf('h''''_rms (RMS curvature in mm^-1): %.4f mm^-1\n', h_rms_double_prime);

% Define a maximum frequency for the cutoff (example value, adjust based on surface size)
q_max = max(q(:)) * 0.9;  % This removes very high-frequency components

% Apply the cutoff to the PSD
psd2_z_cutoff = psd2_z(q <= q_max);

% Recalculate the RMS values with the cutoff applied
h_rms_cutoff = sqrt(sum(psd2_z_cutoff(:)));  % Summing over the allowed frequencies
h_rms_prime_cutoff = sqrt(sum(psd2_z_cutoff(:)));  % Dimensionless
h_rms_double_prime_cutoff = sqrt(sum((q(q <= q_max).^2) .* psd2_z_cutoff(:)));  % Units of mm^-1

% Display the results with the cutoff
fprintf('h_rms (RMS height with cutoff in mm): %.4f mm\n', h_rms_cutoff);
fprintf('h''_rms (RMS slope with cutoff dimensionless): %.4f\n', h_rms_prime_cutoff);
fprintf('h''''_rms (RMS curvature with cutoff in mm^-1): %.4f mm^-1\n', h_rms_double_prime_cutoff);

%% Multiscale Analysis: Power Spectrum at Different Scales
clc
close all

% profile 1
multiscaleAnalysis(xxxq_new1, yyyq_new1, zzzq_new1);

% profile 2
% multiscaleAnalysis(xxxq_new2, yyyq_new2, zzzq_new2);

%% Anisotropy Analysis
clc
close all

% profile 1
analyzeIsotropy(xxxq_new1, yyyq_new1, zzzq_new1);

% profile 2
analyzeIsotropy(xxxq_new2, yyyq_new2, zzzq_new2);

legend('upper profile','lower profile')
xlabel('$$x/x_{max}$$','Interpreter','latex','FontSize',16);
ylabel('$$H$$','Interpreter','latex','FontSize',16)
set(gca,'FontSize',16)
ylim([0 0.7])

%% PSD, C_iso, C_aniso
clc
close all

% profile 1
plotPSD2D_3D(xxxq_new1, yyyq_new1, zzzq_new1);

% profile 2
% plotPSD2D_3D(xxxq_new2, yyyq_new2, zzzq_new2);

%% End