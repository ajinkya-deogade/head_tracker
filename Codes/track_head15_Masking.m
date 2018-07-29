%% Start with nothing (an empty workspace)
clear all;
clc;
% close all;

%% Input/Output Folders
dataFolder = 'C:\Users\adeogade\Desktop\Frames\';
dataDir = strcat(dataFolder, 'Frames_Cropped1\');
filesArray = dir(strcat(dataDir,'*.tiff'));
maxFile = 3000;
cd(dataDir);

%% Load Background Image
bk = strcat(dataDir,'BK\bk.tiff');
ibg = imread(bk);

%% Initialization
fileCount = 0;
flag = 1;
iniFile = 1;
finalFile = maxFile - 50;
cgTime = 1;
skip_frame = 0;
disk_filter = fspecial('disk',2);
se = strel('disk', 2);
boxsize = 10;
boxsize2 = 5;
nextbboxref = nan(finalFile - iniFile + 1, 4);
nextbboxref_head = nan(finalFile - iniFile, 4);
nextbboxref_tail = nan(finalFile - iniFile, 4);
bboxref = nan(finalFile - iniFile + 1, 4);
larvae = struct('centroid','area','perimeter','orientation','majoraxislength','minoraxislength');
larvae.centroid = struct('x','y');
head_number = 0;
tail_number = 0;
head_not_found = 0;
tail_not_found = 0;

for file_number = iniFile:finalFile
    fileCount = fileCount + 1;
    filename = filesArray(file_number).name;
    
    i0 = imread(filename);
    i00 = i0;
    i1 = ibg - i0;
    i2 = imfilter(i1, disk_filter);
    i3 = imadjust(i2);
    if file_number == iniFile
        i4 = im2bw(i3);
    end
    head_found = 0;
    tail_found = 0;
    
    if file_number == iniFile
        larvae_stats_ini = regionprops(i4,'BoundingBox','Area');
        biggest_object_ini = find([larvae_stats_ini.Area] == max([larvae_stats_ini.Area]));
        biggest_box_ini = larvae_stats_ini(biggest_object_ini).BoundingBox;
        point_x = biggest_box_ini(1)-boxsize;
        point_y = biggest_box_ini(2)-boxsize;
        width = biggest_box_ini(3)+2*boxsize;
        height = biggest_box_ini(4)+2*boxsize;
        if point_x < 0, point_x = 0; end
        if point_y < 0, point_y = 0; end
        bbox = [point_x point_y width height];
        nextbboxref(iniFile,:) = bbox;
        nextbboxref(iniFile+1,:) = bbox;
    end
    
    i_crop = imcrop(i3, nextbboxref(file_number,:));
    i_crop = im2bw(i_crop,0.25);
    
    i5 = imclearborder(i_crop);
    i6 = bwareaopen(i5,10);
    i7 = imfill(i6,'holes');
    i8 = bwmorph(i7,'thin', Inf);
    i9 = bwmorph(i8,'spur', 5);
    i10 = bwareaopen(i9, 10);
    
    
    %% Separate out Head and Tail and get the bounding box
    
    i_crop3 = imcrop(i0, nextbboxref(file_number,:));
    i_crop31 = imadjust(i_crop3,[0 0.5], [0 1]);
    i_crop32 = imclose(i_crop3, se);
    %     figure(2), subplot(1,5,2),imshow(i_crop32);
    i_crop35 = i_crop32 - i_crop3;
    i_crop35 = imadjust(i_crop35,[0 0.8],[0 1]);
    %     figure(2), subplot(1,5,3), imshow(i_crop35);
    i_crop36 = im2bw(i_crop35, 0.2);
    %     figure(2), subplot(1,5,4), imshow(i_crop36);
    i_crop37 = imclearborder(i_crop36);
    i_crop37 = bwareaopen(i_crop37, 7, 4);
    i_crop37 = imfill(i_crop37,'holes');
    i_crop37 = bwmorph(i_crop37,'bridge',2);
    %     figure(2), subplot(1,5,5),imshow(i_crop37);
    
    head_stats = regionprops(i_crop37,i_crop31,'Centroid','Area','BoundingBox','Perimeter','Orientation','MajorAxisLength','MinorAxisLength','MeanIntensity','WeightedCentroid','MaxIntensity','MinIntensity');
    area_head = sort([head_stats.Area],'descend');
    intensity_head = sort([head_stats.MeanIntensity],'ascend');
    intensity_min_head = sort([head_stats.MinIntensity],'ascend');
    if (head_number == 0 | tail_number == 0) & (~isempty(area_head))
    head1 = find([head_stats.Area] == area_head(1) & [head_stats.MinIntensity] < 70 )
    elseif head_number > 0 & tail_number > 0
        head1 = find([head_stats.Area] == area_head(1) & [head_stats.MinIntensity] < 70 )
    end
    
    if head_number > 0 & ~isempty(area_head)
        biggest_object2 = find([head_stats.Area] == area_head(1));
        lowest_intensity = find([head_stats.MeanIntensity] == intensity_head(1));
        if biggest_object2 > 0 & lowest_intensity < 110
            head_number = head_number + 1;
            head_found = 1;
            headIntensity(head_number) = head_stats(biggest_object2).MeanIntensity;
            headCentroid(head_number,:) = head_stats(biggest_object2).WeightedCentroid;
            biggest_box2 = head_stats(biggest_object2).BoundingBox;
            bbox_head = biggest_box2;
            point_x_h = biggest_box2(1)-boxsize2;
            point_y_h = biggest_box2(2)-boxsize2;
            width_h = biggest_box2(3)+2*boxsize2;
            height_h = biggest_box2(4)+2*boxsize2;
            if point_x_h < 0, point_x_h = 0; end
            if point_y_h < 0, point_y_h = 0; end
            absolute_bbox_head = [point_x_h point_y_h width_h height_h];
            nextbboxref_head((head_number),:) = absolute_bbox_head;
            dist_Headbox(file_number) = sqrt(abs((nextbboxref_head(head_number,1)-nextbboxref_head(head_number-1,1))^2 - (nextbboxref_head(head_number,2)-nextbboxref_head(head_number-1,2))^2));
        else
        nextbboxref_head((head_number),:) = nextbboxref_head((head_number-1),:);
        headCentroid(head_number,:) = headCentroid(head_number-1,:);
        end
    end
    
    if tail_number > 0 & length(area_head) > 1
        biggest_object3 = find([head_stats.Area] == area_head(2));
        if biggest_object3 > 0
            tail_number = tail_number + 1;
            tail_found = 1;
            tailIntensity(tail_number) = head_stats(biggest_object3).MeanIntensity;
            tailCentroid(tail_number,:) = head_stats(biggest_object3).WeightedCentroid;
            biggest_box3 = head_stats(biggest_object3).BoundingBox;
            bbox_tail = biggest_box3;
            point_x_t = biggest_box3(1)-boxsize2;
            point_y_t = biggest_box3(2)-boxsize2;
            width_t = biggest_box3(3)+2*boxsize2;
            height_t = biggest_box3(4)+2*boxsize2;
            if point_x_t < 0, point_x_t = 0; end
            if point_y_t < 0, point_y_t = 0; end
            absolute_bbox_tail = [point_x_t point_y_t width_t height_t];
            nextbboxref_tail((tail_number),:) = absolute_bbox_tail;
            dist_Tailbox(file_number) = sqrt(abs((nextbboxref_tail(tail_number,1)-nextbboxref_tail(tail_number-1,1))^2 - (nextbboxref_tail(tail_number,2)-nextbboxref_tail(tail_number-1,2))^2));
            if dist_Tailbox(file_number) > 4
                nextbboxref_tail((tail_number),:) = nextbboxref_tail((tail_number-1),:);
            end
        end
    end
    
    if (head_number == 0 | tail_number == 0) & (~isempty(area_head))
        if  head_number == 0
            head_number = head_number + 1;
            head_found = 1;
            biggest_object2 = find([head_stats.Area] == area_head(1));
            biggest_box2 = head_stats(biggest_object).BoundingBox;
            headArea(head_number) = head_stats(biggest_object2).Area;
            headCentroid(head_number,:) = head_stats(biggest_object2).WeightedCentroid;
            point_x_h = biggest_box2(1)-boxsize2;
            point_y_h = biggest_box2(2)-boxsize2;
            width_h = biggest_box2(3)+2*boxsize2;
            height_h = biggest_box2(4)+2*boxsize2;
            if point_x_h < 0, point_x_h = 0; end
            if point_y_h < 0, point_y_h = 0; end
            bbox_head = [point_x_h point_y_h width_h height_h];
            nextbboxref_head(head_number,:) = bbox_head;
            nextbboxref_head(head_number+1,:) = bbox_head;
        end
        
        if  tail_number == 0 & length(area_head) > 1
            tail_number = tail_number + 1;
            tail_found = 1;
            
            biggest_object3 = find([head_stats.Area] == area_head(2));
            biggest_box3 = head_stats(biggest_object3).BoundingBox;
            tailArea(tail_number) = head_stats(biggest_object3).Area;
            tailCentroid(tail_number,:) = head_stats(biggest_object3).WeightedCentroid;
            point_x_t = biggest_box3(1)-boxsize2;
            point_y_t = biggest_box3(2)-boxsize2;
            width_t = biggest_box3(3)+2*boxsize2;
            height_t = biggest_box3(4)+2*boxsize2;
            if point_x_t < 0, point_x_t = 0; end
            if point_y_t < 0, point_y_t = 0; end
            bbox_tail = [point_x_t point_y_t width_t height_t];
            nextbboxref_tail(tail_number,:) = bbox_tail;
            nextbboxref_tail(tail_number+1,:) = bbox_tail;
        end
    end
    
    if head_found == 0 & head_number > 1
        head_not_found = head_not_found + 1;
        head_number = head_number + 1;
        nextbboxref_head(head_number,:) = nextbboxref_head(head_number-1,:);
        headCentroid(head_number,:) = headCentroid(head_number-1,:);
    end
    
    if tail_found == 0 & tail_number > 1
        tail_not_found = tail_not_found + 1;
        tail_number = tail_number + 1;
        nextbboxref_tail(tail_number,:) = nextbboxref_tail(tail_number-1,:);
        tailCentroid(tail_number,:) = tailCentroid(tail_number-1,:);
    end
    
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
    
    if(length(skeleton.endpoints) == 2 & aspectRatio > 1.5)
        skeleton_backbone = skeleton.endpoints(1);
        %         skeleton_backboneLength = length(skeleton_backbone);
        counter = 0;
        while (length(skeleton_backbone) < skeleton.length & counter < 1000)
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
        for counter = 2:skeleton.backboneLength
            s(counter) = s(counter-1) + sqrt((x(counter) - x(counter-1))^2 + (y(counter) - y(counter-1))^2);
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
    end
    %% Head and Tail coordinates in whole frame
    xh_frame(file_number,:) = xh + nextbboxref(file_number,1); yh_frame(file_number,:) = yh + nextbboxref(file_number,2);
    xt_frame(file_number,:) = xt + nextbboxref(file_number,1); yt_frame(file_number,:) = yt + nextbboxref(file_number,2);
    xc_frame(file_number,:) = xc + nextbboxref(file_number,1); yc_frame(file_number,:) = yc + nextbboxref(file_number,2);
    
    %%  Head and Tail Bounding Box in whole frame
    
    if  head_number > 0
        nextbboxref_head_frame(head_number,1) = nextbboxref_head(head_number,1) + nextbboxref(file_number,1);
        nextbboxref_head_frame(head_number,2) = nextbboxref_head(head_number,2) + nextbboxref(file_number,2);
        nextbboxref_head_frame(head_number,3) = nextbboxref_head(head_number,3);
        nextbboxref_head_frame(head_number,4) = nextbboxref_head(head_number,4);
        
        headCentroid_frame(head_number,1) = headCentroid(head_number,1) + nextbboxref(file_number,1);
        headCentroid_frame(head_number,2) = headCentroid(head_number,2) + nextbboxref(file_number,2);
    end
    
    if tail_number > 0
        nextbboxref_tail_frame(tail_number,1) = nextbboxref_tail(tail_number,1) + nextbboxref(file_number,1);
        nextbboxref_tail_frame(tail_number,2) = nextbboxref_tail(tail_number,2) + nextbboxref(file_number,2);
        nextbboxref_tail_frame(tail_number,3) = nextbboxref_tail(tail_number,3);
        nextbboxref_tail_frame(tail_number,4) = nextbboxref_tail(tail_number,4);
        
        tailCentroid_frame(tail_number,1) = tailCentroid(tail_number,1) + nextbboxref(file_number,1);
        tailCentroid_frame(tail_number,2) = tailCentroid(tail_number,2) + nextbboxref(file_number,2);
    end
    
    %% And save current values for subsequent comparison in next frame
    oxh = xh; oyh = yh;
    oxt = xt; oyt = yt;
    oxc = xc; oyc = yc;
    
    %% Plotting the data
    figure(1), subplot(1,4,1), imshow(i0),
    hold on,
    rectangle('Position', nextbboxref(file_number,:),'EdgeColor','w'),
    plot(xh_frame(file_number,1),yh_frame(file_number,1),'r*'),
    plot(xt_frame(file_number,1),yt_frame(file_number,1),'g+'),
    plot(xc_frame(file_number,1),yc_frame(file_number,1),'bd'),
    subplot(1,4,3), imshow(i0),
    hold on,
    if head_number > 0, rectangle('Position', nextbboxref_head_frame(head_number,:),'EdgeColor','r'); plot(headCentroid_frame(head_number,1),headCentroid_frame(head_number,2),'r*'); end
    if tail_number > 0, rectangle('Position', nextbboxref_tail_frame(tail_number,:),'EdgeColor','g'); plot(tailCentroid_frame(tail_number,1),tailCentroid_frame(tail_number,2),'g+'); end
    hold off;
    subplot(1,4,4),imshow(i_crop37),
    hold on,
    if head_number > 0, plot(headCentroid(head_number,1),headCentroid(head_number,2),'r*'); end
    if tail_number > 0, plot(tailCentroid(tail_number,1),tailCentroid(tail_number,2),'g*'); end
    hold off;
end