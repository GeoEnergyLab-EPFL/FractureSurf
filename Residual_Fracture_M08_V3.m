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
image_directory = 'C:\Users\mtale\OneDrive\Desktop\CT-Scan\Sandstone CT-Scan\ID01-240125-M08\SlicesY-ID01-240125-M08';
image_files = dir(fullfile(image_directory, '*.tif')); 

%% Main loop
% Number of figures to include in the path
num_figures = numel(image_files);

Opening_Profile_1=[];

R_m1=[];
R_m2=[];
R_v=[];
DF_m1=[];
DF_m2=[];
DF_v=[];

x_surface = [];         
y_surface1 = []; 
z_surface1 = [];        
y_surface2 = []; 
z_surface2 = [];

% Loop through each image
for i = 2400:1:num_figures-80 % It has been eliminated with the CT-scan data
    
    disp(['i = ' num2str(i) ', ' num2str(i/numel(image_files)*100) ' % of images are processed.']);
    
    % Read the TIFF image
    image = imread(fullfile(image_directory, image_files(i).name));

    % Crop Image
    image_crop=image(130:end-55,:);

    % Normalizing the image
    normalized_image = double(image_crop) / double(intmax('uint16'));

    % Enhancing by Applying a median filter and remove noise
    filteredImage = medfilt2(normalized_image, [3, 3]);
    enhancedImage = adapthisteq(filteredImage,'NumTiles',[8 8],'ClipLimit',0.0005);
    enhancedImage(enhancedImage(:,:)>0.3)=1;
    enhancedImage=(enhancedImage-1)*-1;
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
    for kk=1:size(col_indices,1)
        ind=find(col_indices==col_indices(kk));
        if numel(ind)>1 
            if abs(row_indices(max(ind))-row_indices(min(ind)))<20
                lin1x(zz)=col_indices(kk)*10;
                lin1y(zz)=row_indices(min(ind))*10;
                lin2x(zz)=col_indices(kk)*10;
                lin2y(zz)=row_indices(max(ind))*10;
                midlinx(zz)=col_indices(kk)*10;
                midliny(zz)=(lin1y(zz)+lin2y(zz))/2;
                zz=zz+1;
            end
        end
    end
    
    % Store the x, y, and z values for this loop
    x_surface = [x_surface, lin1x/1000];
    % 1st surface
    y_surface1 = [y_surface1, i * 10 * ones(1, length(lin1x))/1000]; 
    z_surface1 = [z_surface1, lin1y/1000];
    
    % 2nd surface         
    y_surface2 = [y_surface2, i * 10 * ones(1, length(lin2x))/1000]; 
    z_surface2 = [z_surface2, lin2y/1000];
    
    % Opening perpendicular to the mid-surface
    W_pr=NormDistFunc([midlinx',midliny'],[lin1x',lin1y'])+...
         NormDistFunc([midlinx',midliny'],[lin2x',lin2y']);

    W_pr_nonan=W_pr(~isnan(W_pr));
    
    figure(2)
    subplot(2,1,1)
    plot(lin1x,W_pr,'k.','markersize',3);
    yline(median(W_pr_nonan),'m--')
    xlim([0 14500])
    title('fracture opening');
    hold off;
    
    subplot(2,1,2)
    plot(lin1x,lin1y,'r.','markersize',3)
    hold on
    plot(lin2x,lin2y,'b.','markersize',3)
    plot(midlinx,midliny,'k.','markersize',3)
    xlim([0 14500])
    title('detected fracture surface')
    hold off;
    pause(0.1)
    
    Opening_Profile_1=[Opening_Profile_1;[i*10,median(W_pr_nonan)]];
    
    % Roughness Analysis
    [Ra_v, Rq_v, Rz_v] = surfaceRoughnessAnalysis(midliny);
    R_v=[R_v;[i*10,[Ra_v, Rq_v, Rz_v]]];
    
    D_m1=NormDistFunc([midlinx',midliny'],[lin1x',lin1y']);
    D_m1=D_m1(~isnan(D_m1));
    [Ra_m1, Rq_m1, Rz_m1] = surfaceRoughnessAnalysis(D_m1);
    R_m1=[R_m1;[i*10,[Ra_m1, Rq_m1, Rz_m1]]];


    D_m2=NormDistFunc([midlinx',midliny'],[lin2x',lin2y']);
    D_m2=D_m2(~isnan(D_m2));
    [Ra_m2, Rq_m2, Rz_m2] = surfaceRoughnessAnalysis(D_m2);

    R_m2=[R_m2;[i*10,[Ra_m2, Rq_m2, Rz_m2]]];
    

    D = fractalDimension(D_m1);
    DF_m1 =[DF_m1;[i*10,D]];

    D = fractalDimension(D_m2);
    DF_m2 =[DF_m2;[i*10,D]];

    D = fractalDimension(midliny);
    DF_v =[DF_v;[i*10,D]];

end

%% Clean the data - if needed for ploting
% index=find(z_surface1>1.5 & z_surface1<4.3);
% x_surface=x_surface(index);
% y_surface1=y_surface1(index);
% z_surface1=z_surface1(index);
% y_surface2=y_surface2(index);
% z_surface2=z_surface2(index);

%% 3D fracture suface
close all
% 1st surface - cloud points
figure
scatter3(x_surface, (y_surface1-y_surface1(end))*(-1), z_surface1, 10, z_surface1, 'filled');
xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
title('upper surface');
% colorbar;
ax = gca; 
ax.DataAspectRatio = [0.8 0.5 0.3]; 

% 1st surface - surface
xxx = x_surface;
yyy = (y_surface1-y_surface1(end))*(-1);
zzz = z_surface1;
[xxxq, yyyq] = meshgrid(min(xxx):0.1:max(xxx), min(yyy):0.1:max(yyy));
zzzq = griddata(xxx, yyy, zzz, xxxq, yyyq, 'cubic');

figure
surf(xxxq, yyyq, zzzq);
xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
title('upper surface');
ax = gca; 
ax.DataAspectRatio = [0.8 0.5 0.3]; 


% 2st surface - cloud points
% figure
% scatter3(x_surface, (y_surface2-y_surface2(end))*(-1), z_surface2, 10, z_surface2, 'filled');
% xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
% ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
% zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
% title('upper surface');
% % colorbar;
% ax = gca; 
% ax.DataAspectRatio = [0.8 0.5 0.3]; 
% 
% % 2st surface - surface
% xxx = x_surface;
% yyy = (y_surface2-y_surface2(end))*(-1);
% zzz = z_surface2;
% [xxxq, yyyq] = meshgrid(min(xxx):0.1:max(xxx), min(yyy):0.1:max(yyy));
% zzzq = griddata(xxx, yyy, zzz, xxxq, yyyq, 'cubic');
% 
% figure
% surf(xxxq, yyyq, zzzq);
% xlabel('$$x\;(mm)$$','Interpreter','latex','FontSize',16);
% ylabel('$$r\;(mm)$$','Interpreter','latex','FontSize',16);
% zlabel('$$z\;(mm)$$','Interpreter','latex','FontSize',16);
% title('upper surface');
% ax = gca; 
% ax.DataAspectRatio = [0.8 0.5 0.3]; 

%% Residual opening
close all

Opening_Profile_Perpendicular=[((Opening_Profile_1(:,1)-Opening_Profile_1(end,1))*-1)/(10^(4)),Opening_Profile_1(:,2)];
% Opening_Profile_Perpendicular(Opening_Profile_Perpendicular(:,1)>3.05,:)=[];

YY=smooth(Opening_Profile_Perpendicular(:,1),Opening_Profile_Perpendicular(:,2),0.5);

figure
plot(Opening_Profile_Perpendicular(:,1),Opening_Profile_Perpendicular(:,2),'linewidth',1.5,'color',[160 160 160]/255);
hold on;
plot(Opening_Profile_Perpendicular(:,1),YY,'linewidth',2);
xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
ylabel('$$w\;(\mu m)$$','Interpreter','latex','FontSize',16)
set(gca,'FontSize',14)
% xlim([0 3])

figure
plot(Opening_Profile_Perpendicular(:,1),YY,'linewidth',2);
xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
ylabel('$$w\;(\mu m)$$','Interpreter','latex','FontSize',16)
set(gca,'FontSize',14)
% xlim([0 3])

Data_Pro=[Opening_Profile_Perpendicular(:,1), YY];
% Data_Pro(Data_Pro(:,1)>3,:)=[];

error=10*ones(1,size(Data_Pro,1));
figure
e=errorbar(Data_Pro(:,1), Data_Pro(:,2), error, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
e.Color = [220 220 220]/255;
hold on
plot(Data_Pro(:,1), Data_Pro(:,2),'linewidth',1.5,'Color','k');
xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
ylabel('$$w\;(\mu m)$$','Interpreter','latex','FontSize',16)
set(gca,'FontSize',16)
% grid on
% xlim([0 4])
% ylim([20 80])

%% Roughness plots
close all

figure
plot((-1)*(R_v(:,1)-R_v(end,1))/(10^(4)),R_v(:,2))
hold on
plot((-1)*(R_v(:,1)-R_v(end,1))/(10^(4)),R_v(:,3))
plot((-1)*(R_v(:,1)-R_v(end,1))/(10^(4)),R_v(:,4))
legend('Average Roughness','Root Mean Square Roughness','Peak-to-Valley Height')
ylabel('$$\alpha\;(\mu m)$$','Interpreter','latex','FontSize',16)
xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
title('mid-surface')
set(gca,'FontSize',14)
% xlim([0 3])

figure
plot((-1)*(R_m1(:,1)-R_m1(end,1))/(10^(4)),R_m1(:,2))
hold on
plot((-1)*(R_m1(:,1)-R_m1(end,1))/(10^(4)),R_m1(:,3))
plot((-1)*(R_m1(:,1)-R_m1(end,1))/(10^(4)),R_m1(:,4))
legend('Average Roughness','Root Mean Square Roughness','Peak-to-Valley Height')
xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
ylabel('$$\alpha\;(\mu m)$$','Interpreter','latex','FontSize',16)
title('upper surface')
set(gca,'FontSize',14)
% xlim([0 3])

figure
plot((-1)*(R_m2(:,1)-R_m2(end,1))/(10^(4)),R_m2(:,2))
hold on
plot((-1)*(R_m2(:,1)-R_m2(end,1))/(10^(4)),R_m2(:,3))
plot((-1)*(R_m2(:,1)-R_m2(end,1))/(10^(4)),R_m2(:,4))
legend('Average Roughness','Root Mean Square Roughness','Peak-to-Valley Height')
xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
ylabel('$$\alpha\;(\mu m)$$','Interpreter','latex','FontSize',16)
title('lower surface')
set(gca,'FontSize',14)
% xlim([0 3])

figure
plot((-1)*(DF_m1(:,1)-DF_m1(end,1))/(10^(4)),DF_m1(:,2))
xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
ylabel('$$D\;(-)$$','Interpreter','latex','FontSize',16)
set(gca,'FontSize',14)
title('D upper surface')
% xlim([0 3])

figure
plot((-1)*(DF_m2(:,1)-DF_m2(end,1))/(10^(4)),DF_m2(:,2))
xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
ylabel('$$D\;(-)$$','Interpreter','latex','FontSize',16)
set(gca,'FontSize',14)
title('D lower surface')
% xlim([0 3])

figure
plot((-1)*(DF_v(:,1)-DF_v(end,1))/(10^(4)),DF_v(:,2))
xlabel('$$r\;(cm)$$','Interpreter','latex','FontSize',16)
ylabel('$$D\;(-)$$','Interpreter','latex','FontSize',16)
set(gca,'FontSize',14)
title('D mid surface')
% xlim([0 3])

%% End