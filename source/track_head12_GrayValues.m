%% Start with nothing (an empty workspace)
clear all;
% close all;

%% Input/Output Folders
dataFolder = 'C:\Users\adeogade\Desktop\Frames\';
dataDir = strcat(dataFolder, 'Frames_Cropped1\');
filesArray = dir(strcat(dataDir,'*.tiff'));
maxFile = 3000;
cd(dataDir);

%% Load Background Image
bk = strcat(dataDir, 'BK\bk.tiff');
ibg = imread(bk);

%% Initialization
fileCount = 0;
flag = 1;
iniFile = 1;
finalFile = maxFile - 50;
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
    %     i0 = rgb2gray(i0);
    i00 = i0;
    i1 = ibg - i0;
    i2 = imfilter(i1, disk_filter);
    i3 = imadjust(i2);
    if file_number == iniFile
    i4 = im2bw(i3);
    end
    
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
    
    
    i_crop = imcrop(i3, nextbboxref(file_number,:));
    if file_number > iniFile
        i_crop = im2bw(i_crop,0.25);
    end
    
    i_crop2 = imcrop(i0, nextbboxref(file_number,:));
    
    
    i5 = imclearborder(i_crop);
    i6 = bwareaopen(i5,10);
    i7 = imfill(i6,'holes');
    i8 = bwmorph(i7,'thin', Inf);
    i9 = bwmorph(i8,'spur', 5);
    i10 = bwareaopen(i9, 10);
    
%     das = im2bw(i_crop2,0.2);
    
    figure (2),
    subplot(1,3,1), imshow(i_crop2);
    subplot(1,3,2), imshow(i5);
    
    
    %     subplot(2,2,2), imshow(i6);
    %     subplot(2,2,3), imshow(i8);
    %     subplot(2,2,4), imshow(i10);
    
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
    
    [rMain, cMain] = find(bwlabel(i7) == biggest_object);
    iL = zeros(size(i7,1),size(i7,2));
    for k=1:length(rMain)
        iL(rMain(k),cMain(k))=1;
    end
    iL = imfill(iL,'holes');
    iperim = bwperim(iL);
    centerIndp = find(iperim>0);
    [py,px] = ind2sub(size(iperim),centerIndp);
    contour = bwtraceboundary(iperim, [py(5) px(5)],'N',8, Inf,'counterclockwise');
    roi = roipoly(i_crop2, contour(:,2), contour(:,1));
%     H = fspecial('unsharp');
%     J = roifilt2(H, roi, i_crop2);
 cop = i_crop2;
 cop(~roi) = 0;
    figure(2), subplot(1,3,3), imshow(cop);
    
%     [Ys, Xs] = smooth_contours(contour(:,2), contour(:,1), 10);
%     contour(:,1) = Xs;
% %     contour(:,2) = Ys;
% for i = 1:length(contour)
%     i_crop2(contour(i,1),contour(i,2)) = 1;
% end
%         subplot(1,3,3), imshow(i_crop2);
    %% Skeleton Information of Larvae
    
    skeleton.index = find(i10 > 0);
    skeleton.length = length(skeleton.index);
    [skeleton_y, skeleton_x] = ind2sub(size(i10),skeleton.index);
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
            skeleton_backbone = cat(2, skeleton_backbone, max(next_index));
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
        
        %% The Proximity Rule
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
        
        %% Head and Tail coordinates in whole frame
        xh_frame(file_number,:) = xh + nextbboxref(file_number,1); yh_frame(file_number,:) = yh + nextbboxref(file_number,2);
        xt_frame(file_number,:) = xt + nextbboxref(file_number,1); yt_frame(file_number,:) = yt + nextbboxref(file_number,2);
        xc_frame(file_number,:) = xc + nextbboxref(file_number,1); yc_frame(file_number,:) = yc + nextbboxref(file_number,2);
        
        %% Define Head and Tail Region
        [x0, y0, coordinate_real] = signature(contour, xc, yc);
        st_select = find(coordinate_real(:,6) < 25);
        h = zeros(size(iperim));
        for k = 1:length(st_select)
            h(coordinate_real(st_select(k),1),coordinate_real(st_select(k),2)) = 1;
        end
        h = bwmorph(h, 'bridge');
        %         figure(1),
        %         subplot(1,3,1), imshow(i5);
        %         subplot(1,3,2), imshow(h);
        %         subplot(1,3,3), plot(coordinate_real(st_select,3), coordinate_real(st_select,4),'r.');
        
        %%   Directional Information of Movement
        new_head_orient(file_number,1) = atan2(yh_frame(file_number,1),xh_frame(file_number,1));
        new_head_mag(file_number,1) = sqrt(yh_frame(file_number,1)^2 + xh_frame(file_number,1)^2);
        new_head_magx = new_head_mag(file_number,1)*cos(new_head_orient(file_number,1));
        new_head_magy = new_head_mag(file_number,1)*sin(new_head_orient(file_number,1));
        if file_number > iniFile
            old_head_magx = new_head_mag(file_number-1,1)*cos(new_head_orient(file_number-1,1));
            old_head_magy = new_head_mag(file_number-1,1)*sin(new_head_orient(file_number-1,1));
            directionChangex = new_head_magx + old_head_magx;
            directionChangey = new_head_magy + old_head_magy;
            directionChange_mag(file_number,1) = sqrt(directionChangex^2 + directionChangey^2);
            directionChange_orient(file_number,1) = atan2(directionChangey,directionChangex);
        end
        
        %         if file_number > iniFile+1
        %             figure(1), subplot(1,3,1), plot((iniFile+1:file_number),directionChange_mag((iniFile+1:file_number),1)), title('Change in Magnitude');
        %             subplot(1,3,2), plot((iniFile+1:file_number),rad2deg(directionChange_orient((iniFile+1:file_number),1))), title('Change in Direction');
        %             subplot(1,3,3), imshow(i5);
        %         end
        %%  Define the coarse directions
        %         if
        
        %% Cut the neck according to the direction of heading
        
        larvaeH = i7;
        %                 larvaeH(yh, xh) = 0;
        
        
        %         k = 1;
        %         for i=1:5
        %         point_right(i) =
        %         while (larvaeH(yh, xh + k) == 255 || larvaeH(yh, xh - k) == 255) && xh+k <= length(larvaeH(yh,:))
        %             larvaeH(yh, xh + k) = 0;
        %             larvaeH(yh, xh - k) = 0;
        %             k = k+1;
        %         end
        %         figure(1), imshow(larvaeH);
        %% Get the grey values in Head region
        
        
        
        %% And save current values for subsequent comparison in next frame
        oxh = xh; oyh = yh;
        oxt = xt; oyt = yt;
        oxc = xc; oyc = yc;
        
        %% Calculating the distance of head from the centre of body
        distance_head_center(file_number,1) = sqrt((xh_frame(file_number,:) - xc_frame(file_number,:))^2 + (yh_frame(file_number,:)-yc_frame(file_number,:))^2);
        
        %% Calculating the distance of head from previous head position
        if file_number > iniFile
            distance_head_previous(file_number,1) = sqrt((xh_frame(file_number,:) - xh_frame(file_number-1,:))^2 + (yh_frame(file_number,:)-yh_frame(file_number-1,:))^2);
            if distance_head_previous > 20  continue; end
        end
        
        %         %% Plotting the data
        %                 figure(3), subplot(1,3,1), imshow(i4),
        %                 hold on,
        %                 rectangle('Position', nextbboxref(file_number,:),'EdgeColor','w'),
        %                 plot(xh_frame(file_number,1),yh_frame(file_number,1),'r*'),
        %                 plot(xt_frame(file_number,1),yt_frame(file_number,1),'g+'),
        %                 plot(xc_frame(file_number,1),yc_frame(file_number,1),'bd'),
        %                 hold off;
        %
        %                 figure(3),
        %                 subplot(1,3,2), plot((1:file_number),distance_head_center(:,1)), title('Distance Between Head and Body Center');
        %                 if file_number > 1
        %                     subplot(1,3,3), plot((1:file_number), distance_head_previous(:,1)), title('Distance Between Current Head Position and Previous Head Position');
        %                 end
        
    end
end
sprintf('Total Number of frames skipped : %d', skip_frame);