%% Start with nothing (an empty workspace)
clear all;
close all;

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
iniFile = 1;
finalFile = maxFile;
cgTime = 1;
skip_frame = 0;
Laplacian_of_Gaussian = fspecial('log', [3 3],0.75);
disk_filter = fspecial('disk',2);
boxsize = 10;
head_point = nan(finalFile-iniFile,2);
tail_point = nan(finalFile-iniFile,2);

for file_number = iniFile:finalFile
    
    % A useful counter
    fileCount = fileCount + 1;
    clear x y;
    
    % read each frame from the folder
    filename = filesArray(file_number).name;
    i0 = imread(filename);
    i0 = rgb2gray(i0);
    i00 = i0;
    i3 = ibg - i0;
    %g = graythresh(i1)
    i2 = imfilter(i3, disk_filter);
    i3 = imadjust(i2);
    i4 = im2bw(i3);
    
    bbox_prop = regionprops(i4,'BoundingBox','Area');
    biggest_object = find([bbox_prop.Area] == max([bbox_prop.Area]));
    biggest_box = bbox_prop(biggest_object).BoundingBox;
    bbox = [biggest_box(1)-boxsize biggest_box(2)-boxsize biggest_box(3)+2*boxsize biggest_box(4)+2*boxsize ];
    nextbboxref(1,:) = bbox;
        
    if file_number == 1
        figure(1),imshow(i4), hold on, rectangle('Position', nextbboxref(fileCount,:),'EdgeColor','w'), hold off;
        i_crop = imcrop(i4, nextbboxref(fileCount,:));
%         i_crop = i2;
    else
        i_crop = imcrop(i4, nextbboxref(fileCount,:));
        figure(1),imshow(i4), hold on, rectangle('Position', nextbboxref(fileCount,:),'EdgeColor','w'), hold off;
    end 
    
    i5 = bwmorph(i_crop, 'skel', Inf);
    i5 = bwmorph(i5, 'bridge', 3);
    i5 = bwmorph(i5, 'spur', 5);
    E = bwmorph(i5, 'endpoints');
    [y, x] = find(E);
    head_point(file_number,:) = [y(1,1) x(1,1)];
    tail_point(file_number,:) = [y(2,1) x(2,1)];
%     figure(2), imshow(i5);
%, plot(head_point(file_number,2),head_point(file_number,1),'r*')
           
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
    
    %% Bounding Box for larvae
    boxsize = 10;
    if fileCount == 1
        bbox_prop = regionprops(i4,'BoundingBox','Area');
        biggest_object = find([bbox_prop.Area] == max([bbox_prop.Area]));
        biggest_box = bbox_prop(biggest_object).BoundingBox;
        bbox = biggest_box;
        nextbboxref(1,:) = bbox;
    else
        bbox_prop = regionprops(i4,'BoundingBox','Area');
        biggest_object = find([bbox_prop.Area] == max([bbox_prop.Area]));
        biggest_box = bbox_prop(biggest_object).BoundingBox;
        bbox = biggest_box;
    end
    bboxref(fileCount,:) = bbox;
    absolute_bbox = [bboxref(fileCount,1)-boxsize bboxref(fileCount,2)-boxsize (bboxref(fileCount,3)+2*boxsize) (bboxref(fileCount,4)+2*boxsize)];
    nextbboxref((fileCount+1),:) = absolute_bbox;
    
    %% Find Curvature
    neighbour_distance = 3;
    length_contour = length(contour_body(:,1));
    curVec = nan(1,length_contour);
    for i = 1:length_contour
        index1 = i + neighbour_distance;
        index2 = i - neighbour_distance;
        
        if index1 > length_contour
            index1 = index1 - length_contour;
        end
        if index2 < 1
            index2 = length_contour + index2;
        end
        
        p1(1,1) = contour_body(i,1); p1(1,2) = contour_body(i,2);
        p2(1,1) = contour_body(index1,1); p2(1,2) = contour_body(index1,2);
        p3(1,1) = contour_body(index2,1); p3(1,2) = contour_body(index2,2);
        
        angle1 = atan2((p2(1,1)-p1(1,1)),(p2(1,2) -p1(1,2)))-atan2((p3(1,1)-p1(1,1)),(p3(1,2) -p1(1,2)));
        curVec(i) = angle1;
    end
    curVec = unwrap(curVec);%- mean(unwrap(curVec));
    curVec_degree = rad2deg(curVec);
    circshift(curVec,[0 -neighbour_distance]);
    
    %% Find the maximum positive value: should be the head or tail
    [valMax1, idxMax1] = max(curVec);
    idxMax = idxMax1(1);
    curVec2 = curVec;
    curVec2(idxMax1) = 0;
    length1 = length(contour_body);
    curVec2(length1) = 0;
    
    %% Finally, we find the second maximum
    distCurv = 3;
    length1 = length(contour_body);
    index2 = index1 + distCurv;
    if index2 <= length1
        curVec2(index1:index2) = 0;
    elseif index2 > length1
        indexMove = length1 - index1;
        indexNo = distCurv - indexMove;
        curVec2(index1:length1) = 0;
        curVec2(1:indexNo) = 0;
    end
    
    index3 = index1-distCurv;
    if index3 >= 1
        curVec2(index3:index1) = 0;
    elseif index3 < 1
        indexNo = length1+index3;
        curVec2(1:index1)= 0;
        curVec2(indexNo:length1) = 0;
    end
    [valMax2, idxMax2] = max(curVec2);
    idxMax2 = idxMax2(1);
    
    %% Assign Head and and Tail coordinates
    zhead = [contour_body(idxMax,2) contour_body(idxMax,1)];
    ztail = [contour_body(idxMax2,2) contour_body(idxMax2,1)];
    
    %% Current guesses:
    xt = ztail(2); yt = ztail(1);
    xh = zhead(2); yh = zhead(1);
    
    %% Apply "the proximity rule"
    
    if file_number == iniFile
        xt = ztail(2); yt = ztail(1);
        xh = zhead(2); yh = zhead(1);
    else
        dist_oldHead_newHead = sqrt((xh-oxh)^2 +  (yh-oyh)^2)
        dist_oldTail_newHead = sqrt((xh-oxt)^2 +  (yh-oyt)^2);
        
        
        dist_oldHead_newTail = sqrt((xt-oxh)^2 +  (yt-oyh)^2);
        dist_oldTail_newTail = sqrt((xt-oxt)^2 +  (yt-oyt)^2)
        
        dist_threshold_newHead = dist_oldHead_newHead - dist_oldTail_newHead;
        %         display(dist_threshold_newHead);
        dist_threshold_newTail = dist_oldHead_newTail - dist_oldTail_newTail;
        %         display(dist_threshold_newTail);
        
        if dist_oldHead_newHead < 12
            xh = zhead(2); yh = zhead(1);
        else
            xh = oxh; yh = oyh;
        end
        if dist_oldTail_newTail < 12
            xt = ztail(2); yt = ztail(1);
        else
            xt = oxt; yt = oyt;
        end
        
        % find the distance and swap head and tail if necessary:
        %         if  dist_oldHead_newHead < dist_oldHead_newTail && dist_oldTail_newTail < dist_oldTail_newHead
        %             xh = zhead(2); yh = zhead(1);
        %             xt = ztail(2); yt = ztail(1);
        %         elseif dist_oldHead_newHead > dist_oldHead_newTail && dist_oldTail_newTail < dist_oldTail_newHead
        %             xh = oxh; yh = oyh;
        %             xt = ztail(2); yt = ztail(1);
        %         elseif dist_oldHead_newHead < dist_oldHead_newTail && dist_oldTail_newTail > dist_oldTail_newHead
        %             xh = zhead(2); yh = zhead(1);
        %             xt = oxt; yt = oyt;
        %         elseif dist_oldHead_newHead < dist_oldHead_newTail && dist_oldTail_newTail > dist_oldTail_newHead
        %             xh = ztail(2); yh = ztail(1);
        %             xt = zhead(2); yt = zhead(1);
        %         else
        %             sprintf('Skipped Frame number: %d',file_number)
        %             skip_frame = skip_frame + 1;
        %             oxh = oxh; oyh = oyh;
        %             oxt = oxt; oyt = oyt;
        %           continue
        %         end
    end
    %% And save current values for subsequent comparison in next frame
    oxh = xh; oyh = yh;
    oxt = xt; oyt = yt;
    
    %% Plotting the data
%         if file_number > 1
%             %         ;
%             figure(3), plot(contour_body(:,1), contour_body(:,2), 'r.'), axis equal;
%             hold on;
%             plot(xh, yh,'*g');
%             % plot(xc,yc,'*w');
%             plot(xt, yt,'*y');
%             hold off;
%             %         rectangle('Position',bbox);
%     
%             %         display('Paused: press any key to continue')
%             %         pause
% %             pause;
%             %         axis close;
%     
%         end
end
sprintf('Total Number of frames skipped : %d', skip_frame)

%% Curvature
%     for i = 2:length(contour_body(:,1))
%         if i+1 < length(contour_body(:,1))
%             dxi = contour_body(i,2)-contour_body(i-1,2);
%             dyi = contour_body(i,1)-contour_body(i-1,1);
%             ddxi = contour_body(i-1,2)-2*contour_body(i,2)+contour_body(i+1,2);
%             ddyi = contour_body(i-1,1)-2*contour_body(i,1)+contour_body(i+1,1);
%             curvature(i)=(dxi*ddyi - dyi*ddxi)/((dxi^2 + dyi^2)^1.5);
% %             curvature_degree(i) = rad2deg(curvature(i));
%         else
%             dxi = contour_body(i,2)-contour_body(i-1,2);
%             dyi = contour_body(i,1)-contour_body(i-1,1);
%             ddxi = contour_body(i-1,2)-2*contour_body(i,2)+contour_body(1,2);
%             ddyi = contour_body(i-1,1)-2*contour_body(i,1)+contour_body(1,1);
%             curvature(i)=(dxi*ddyi - dyi*ddxi)/((dxi^2 + dyi^2)^1.5);
% %             curvature_degree(i) = rad2deg(curvature(i));
%         end
%     end
%     curvature = unwrap(curvature) -mean(unwrap(curvature)) ;
%     curvature_degree = unwrap(curvature_degree) - mean(unwrap(curvature));
%     for j = 1:contour_info(2,3)
%         contour_head(:,j) = i_contour(:, contour_info(2,1)+j);
%     end
%     contour_head = contour_head';
%
%     figure(3), plot(contour_head(:,1), contour_head(:,2)), axis equal;
%     hold on;
%     plot(contour_body(:,1), contour_body(:,2));



% max_index = contour_info(1,));


% figure (3);
% plot(i_contour');


%% Fitting
%     contour_body = contour_body';
%     originalSpacing = 1 : length(contour_body(1,:));
%     finerSpacing = 1 : 0.1 : 175;
%     splineXY = spline(originalSpacing, contour_body, finerSpacing);
%     contour_body_fit(1,:) = splineXY(1,:);
%     contour_body_fit(2,:) = splineXY(2,:);
%     contour_body_fit = contour_body_fit';
%     for i = floor((0.95*length(contour_body_fit(1,:)))):length(contour_body_fit(1,:))
%         contour_body_fit(i,1) = contour_body_fit(1,1);
%         contour_body_fit(i,2) = contour_body_fit(1,2);
%     end
