%% Start with nothing (an empty workspace)
clear all;
% close all;

%% Input/Output Folders
dataDir0 = 'F:\Video\';
dirOut2 = 'F:\Video_output';
mkdir(dirOut2);
dataFolder = 'larvae\';
dataDir = strcat(dataDir0,dataFolder,'frames\');
outDir = dataDir;
filesArray = dir(strcat(dataDir,'*.jpeg'));
maxFile = length(filesArray);
cd(dataDir);

%% Load Background Image
ibg = imread('F:\Video\larvae\frames\BK\0000.jpeg');
ibg = ibg(:,:,1);

%% Initialization
fileCount = 0;
flag = 1;
iniFile = 500;
finalFile = maxFile;
cgTime = 1;
skip_frame = 0;
Laplacian_of_Gaussian = fspecial('log', [3 3],0.75);
disk_filter = fspecial('disk',2);
boxsize = 10;
larvae_head = nan(finalFile-iniFile,2);
larvae_tail = nan(finalFile-iniFile,2);
larvae_centriod = nan(finalFile-iniFile,2);
larvae = struct('centroid','area','perimeter','orientation','majoraxislength','minoraxislength');
larvae.centroid = struct('x','y');

for file_number = iniFile:finalFile
    
    % A useful counter
    fileCount = fileCount + 1;
    
    % read each frame from the folder
    filename = filesArray(file_number).name;
    i0 = imread(filename);
    i0 = rgb2gray(i0);
    i00 = i0;
    i1 = ibg - i0;
%     i2 = imfilter(i1, disk_filter);
    i3 = imadjust(i1);
    i4 = im2bw(i3);
    
    if file_number == iniFile
        larvae_stats = regionprops(i4,'BoundingBox','Area');
        biggest_object = find([larvae_stats.Area] == max([larvae_stats.Area]));
        biggest_box = larvae_stats(biggest_object).BoundingBox;
        point_x = biggest_box(1)-boxsize;
        point_y = biggest_box(2)-boxsize;
        width = biggest_box(3)+2*boxsize;
        height = biggest_box(4)+2*boxsize;
        if point_x < 0 point_x = 0; end
        if point_y < 0 point_y = 0; end
        bbox = [point_x point_y width height];
        nextbboxref(iniFile,:) = bbox;
        nextbboxref(iniFile+1,:) = bbox;
    end
    
%     figure(1),imshow(i4), hold on, rectangle('Position', nextbboxref(file_number,:),'EdgeColor','w'), hold off;
    
    i_crop = imcrop(i4, nextbboxref(file_number,:));
    
    i5 = imclearborder(i_crop);
    i6 = bwareaopen(i5,10);
    i7 = imfill(i6,'holes');
    i8 = bwmorph(i7,'thin', Inf);
    i9 = bwmorph(i8,'spur', 6);
%     i10 = bwareaopen(i9, 10);

    figure (2)
    subplot(2,2,1), imshow(i5);
    subplot(2,2,2), imshow(i6);
    subplot(2,2,3), imshow(i8);
    subplot(2,2,4), imshow(i9);
        
    %% Extract Body Contour
    clear contour_info contour_body contour_head contour_body_fit;
    contour_info = []; contour_body = [];
    figure(3);
    i_contour = imcontour(i_crop, 1); axis equal;
    
    contour_start = find(i_contour(1,:) < 1);
    num_of_contours = length(contour_start);
    contour_info = zeros(num_of_contours, 4);
    
    for i = 1:num_of_contours
        contour_info(i,1) = contour_start(i); % Contour Start Points
        if i < num_of_contours
            contour_info(i,2) = contour_start(i+1)-1; % Contour End Points
        else
            contour_info(num_of_contours,2) = num_of_contours;
        end
        contour_info(i,3) = i_contour(2,contour_start(i)); % Contour Length
        contour_info(i,4) = i_contour(1,contour_start(i));
    end
    contour_info';
    contour_info = sortrows(contour_info, -3); % Sort according to rows
    
    for j = 1:contour_info(1,3)
        contour_body(:,j) = i_contour(:, contour_info(1,1) + j); % Pick the biggest contour as larvae
    end
    
    contour_body = contour_body';
    
    %% Complete the Contour
    l = length(contour_body(:,2));
    if contour_body(end,1) ~= contour_body(1,1) || contour_body(end,2) ~= contour_body(1,2)
        for k = 1:20
            contour_body(l+k,1) = contour_body(1,1);
            contour_body(l+k,2) = contour_body(1,2);
        end
        contour_body(l+1,1) = contour_body(1,1);
        contour_body(l+1,2) = contour_body(1,2);
    end
    
    %% Smoothen the Contour
    [Ys, Xs] = smooth_contours(contour_body(:,1), contour_body(:,2), 10);
    contour_body(:,1) = Ys;
    contour_body(:,2) = Xs;
    figure(5), subplot(1,2,1), plot(contour_body(:,2),contour_body(:,1),'g.');
    
    %% Fourier Decomposition
    contour_decompose = fourier_decompose(contour_body, 100, 4);
    
    %% Fourier Reconstruct
    contour_reconstruct = fourier_reconstruct(contour_decompose, 100);
    figure(5), subplot(1,2,2), plot(contour_reconstruct(:,2),contour_reconstruct(:,1),'b.');
    
    %% Morphological Information And Bounding Box for larvae
   
    larvae_stats = regionprops(i_crop,'Centroid','Area','BoundingBox','Perimeter','Orientation','MajorAxisLength','MinorAxisLength');
    biggest_object = find([larvae_stats.Area] == max([larvae_stats.Area]));
    biggest_box = larvae_stats(biggest_object).BoundingBox;
    larvae.centroid.x = larvae_stats(biggest_object).Centroid(1);
    larvae.centroid.y = larvae_stats(biggest_object).Centroid(2);
    larvae.area = larvae_stats(biggest_object).Area;
    larvae.perimeter = larvae_stats(biggest_object).Perimeter;
    larvae.orientation = larvae_stats(biggest_object).Orientation;
    larvae.majoraxislength = larvae_stats(biggest_object).MajorAxisLength;
    larvae.minoraxislength = larvae_stats(biggest_object).MinorAxisLength;
    
    if file_number > iniFile
        bbox = biggest_box;
        bboxref(file_number,:) = bbox;
        point_x = bboxref(file_number,1)+nextbboxref((file_number),1)-boxsize;
        point_y = bboxref(file_number,2)+nextbboxref((file_number),2)-boxsize;
        width = (bboxref(file_number,3)+2*boxsize);
        height = (bboxref(file_number,4)+2*boxsize);
        if point_x < 0 point_x = 0; end
        if point_y < 0 point_y = 0; end
        absolute_bbox = [point_x point_y width height];
        nextbboxref((file_number+1),:) = absolute_bbox;
    end  

    %% Skeleton Information of Larvae
    
    skeleton.index = find(i9 > 0);
    skeleton.length = length(skeleton.index);
    [skeleton_y, skeleton_x] = ind2sub(size(i9),skeleton.index);
    connect = zeros(1,skeleton.length);
    for point = 1:skeleton.length
        i = skeleton_x(point);
        j = skeleton_y(point);
        connect(point) = length(find(skeleton_x >= (i-1) & skeleton_x <= (i+1) & skeleton_y >= (j-1) & skeleton_y <= (j+1)));
    end
    skeleton.endpoints = find(connect == 2);
    length_thresh = 10;
    aspectRatio = larvae.majoraxislength/larvae.minoraxislength;
    
    if(length(skeleton.endpoints) == 2 && aspectRatio > 1.5)
        skeleton_backbone = skeleton.endpoints(1);
        %         skeleton_backboneLength = length(skeleton_backbone);
        counter = 0;
        while (length(skeleton_backbone) < skeleton.length && counter < 1000)
            counter = counter + 1;
            i = skeleton_x(skeleton_backbone(end));
            j = skeleton_y(skeleton_backbone(end));
            vv = find(skeleton_x >= (i-1) & skeleton_x <= (i+1) & skeleton_y >= (j-1) & skeleton_y <= (j+1));
            next_index = setdiff(vv ,skeleton_backbone);
            skeleton_backbone = cat(2, skeleton_backbone, next_index);
        end
        skeleton.backbone = [skeleton_y(skeleton_backbone), skeleton_x(skeleton_backbone)];
        y = skeleton.backbone(:,1); y = reshape(y,1,length(y));
        x = skeleton.backbone(:,2); x = reshape(x,1,length(x));
        skeleton.backboneLength = length(skeleton.backbone);
        s = zeros(1,skeleton.backboneLength);
        for i = 2:skeleton.backboneLength
            s(i) = s(i-1) + sqrt((x(i) - x(i-1))^2 + (y(i) - y(i-1))^2);
        end
        
        [val,idx] = min(abs(s-s(end)/2));
        skeletonMid = skeleton.backbone(idx,:);
        xc = skeletonMid(2);
        yc = skeletonMid(1);
        
        %% Current guesses:
        xt = skeleton_x(skeleton.endpoints(2)); yt = skeleton_y(skeleton.endpoints(2));
        xh = skeleton_x(skeleton.endpoints(1)); yh = skeleton_y(skeleton.endpoints(1));
        
        %% Apply "the proximity rule"
        if file_number > iniFile
            distToTn = sqrt((xt-oxt)^2 + (yt-oyt)^2);
            distToHn = sqrt((xh-oxt)^2 + (yh-oyt)^2);
            endpointHead = skeleton.endpoints(1);
            if distToTn > distToHn
                xh = skeleton_x(skeleton.endpoints(2)); yh = skeleton_y(skeleton.endpoints(2));
                xt = skeleton_x(skeleton.endpoints(1)); yt = skeleton_y(skeleton.endpoints(1));
                endpointHead  = skeleton.endpoints(2);
            end
        end
        
        %% And save current values for subsequent comparison in next frame
        oxh = xh; oyh = yh;
        oxt = xt; oyt = yt;
        
        %% Head and Tail coordinates in whole frame
        xh_frame(file_number,:) = xh + nextbboxref(file_number,1); yh_frame(file_number,:) = yh + nextbboxref(file_number,2);
        xt_frame(file_number,:) = xt + nextbboxref(file_number,1); yt_frame(file_number,:) = yt + nextbboxref(file_number,2);
        xc_frame(file_number,:) = xc + nextbboxref(file_number,1); yc_frame(file_number,:) = yc + nextbboxref(file_number,2);
        
        %% Calculating the distance of head from the centre of body
%         pause
        
        %% Plotting the data
        figure(4),imshow(i4),
        hold on, 
        rectangle('Position', nextbboxref(file_number,:),'EdgeColor','w'), 
        plot(xh_frame(file_number,1),yh_frame(file_number,1),'r*'), 
        plot(xt_frame(file_number,1),yt_frame(file_number,1),'g+'),
        plot(xc_frame(file_number,1),yc_frame(file_number,1),'bd'),
        hold off;
        
    end
end
sprintf('Total Number of frames skipped : %d', skip_frame);