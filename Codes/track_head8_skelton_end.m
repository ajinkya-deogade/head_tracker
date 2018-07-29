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
larvae_head = nan(finalFile-iniFile,2);
larvae_tail = nan(finalFile-iniFile,2);
larvae_centriod = nan(finalFile-iniFile,2);

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
    nextbboxref(2,:) = bbox;
    nextbboxref2(1,:) = bbox;%[biggest_box(1)-50 biggest_box(2)-50 biggest_box(3)+2*50 biggest_box(4)+2*50 ];
    nextbboxref2(2,:) = bbox;
    
%     if file_number == 1
        %         figure(1),imshow(i4), hold on, rectangle('Position', nextbboxref(fileCount,:),'EdgeColor','w'), hold off;
%         i_crop = imcrop(i4, nextbboxref2(file_number,:));
%                 i_crop = i4;
%     else
        i_crop = imcrop(i4, nextbboxref(file_number,:));
        %         figure(1),imshow(i4), hold on, rectangle('Position', nextbboxref(fileCount,:),'EdgeColor','w'), hold off;
%     end
    
    i5 = bwmorph(i_crop, 'skel', Inf);
    i5 = bwmorph(i5, 'bridge', 3);
    i5 = bwmorph(i5, 'spur', 5);
    E = bwmorph(i5, 'endpoints');
    [y, x] = find(E);
    larvae_head(file_number,:) = [y(1,1) x(1,1)];
    larvae_tail(file_number,:) = [y(2,1) x(2,1)];
    %     figure(2), imshow(i5);
    %,
    
    %% Bounding Box for larvae
    boxsize = 10;
%     if file_number == 1
%         bbox_prop = regionprops(i4,'BoundingBox','Area','Centroid');
%         biggest_object = find([bbox_prop.Area] == max([bbox_prop.Area]));
%         biggest_box = bbox_prop(biggest_object).BoundingBox;
%         bbox = biggest_box;
%         nextbboxref(1,:) = bbox;
%     else
        bbox_prop = regionprops(i4,'BoundingBox','Area');
        biggest_object = find([bbox_prop.Area] == max([bbox_prop.Area]));
        biggest_box = bbox_prop(biggest_object).BoundingBox;
        bbox = biggest_box;
%     end
  
    bboxref(file_number,:) = bbox;
    absolute_bbox = [bboxref(file_number,1)+nextbboxref((file_number),1)-boxsize bboxref(file_number,2)+nextbboxref((file_number),2)-boxsize (bboxref(file_number,3)+2*boxsize) (bboxref(file_number,4)+2*boxsize)];
    nextbboxref((file_number+1),:) = absolute_bbox;
    
%     absolutebbox_2 = [bboxref(file_number,1)+nextbboxref2(file_number,1) bboxref(file_number,2)+nextbboxref(file_number,2) bboxref(file_number,3) bboxref(file_number,4)];
%     finalbbox(file_number,:) = absolutebbox_2;
%     nextbboxref2((file_number+1),:) = absolutebbox_2;
    
    %% Current guesses:
    xt = larvae_tail(file_number,2); yt = larvae_tail(file_number,1);
    xh = larvae_head(file_number,2); yh = larvae_head(file_number,1);
    
    %% Apply "the proximity rule"
    if file_number == iniFile
        xt = larvae_tail(1,2); yt = larvae_tail(1,1);
        xh = larvae_head(1,2); yh = larvae_head(1,1);
    else
        distToTn = sqrt((xt-oxt)^2 + (yt-oyt)^2);
        distToHn = sqrt((xt-oxh)^2 + (yt-oyh)^2);
        
        if distToTn > distToHn
            xh = larvae_tail(1,2); yh = larvae_tail(1,1);
            xt = larvae_head(1,2); yt = larvae_head(1,1);
        else
            xh = oxh; yh = oyh;
            xt = oxt; yt = oyt;
            skip_frame = skip_frame + 1;
        end
    end
% if file_number == iniFile
%     xt = larvae_tail(1,2); yt = larvae_tail(1,1);
%     xh = larvae_head(1,2); yh = larvae_head(1,1);
% else
%     dist_oldHead_newHead(file_number) = sqrt((xh-oxh)^2 +  (yh-oyh)^2);
%     dist_oldTail_newHead = sqrt((xh-oxt)^2 +  (yh-oyt)^2);
%     
%     dist_oldHead_newTail = sqrt((xt-oxh)^2 +  (yt-oyh)^2);
%     dist_oldTail_newTail(file_number) = sqrt((xt-oxt)^2 +  (yt-oyt)^2);
    
%     dist_threshold_newHead = dist_oldHead_newHead - dist_oldTail_newHead;
%     dist_threshold_newTail = dist_oldHead_newTail - dist_oldTail_newTail;
%         
%     if dist_oldHead_newHead(file_number) < 12
%         xh = larvae_head(1,2); yh = larvae_head(1,1);
%     else
%         xh = larvae_tail(1,2); yh = larvae_tail(1,1);
%     end
%     if dist_oldTail_newTail(file_number) < 12
%         xt = larvae_tail(1,2); yt = larvae_tail(1,1);
%     else
%         xt = larvae_head(1,2); yt = larvae_head(1,1);
%     end
% end
%% And save current values for subsequent comparison in next frame
oxh = xh; oyh = yh;
oxt = xt; oyt = yt;

%% Head and Taill coordinates in whole frame
xh_frame(file_number,:) = xh + absolute_bbox(1,1); yh_frame(file_number,:) = yh + absolute_bbox(1,2);
xt_frame(file_number,:) = xt + absolute_bbox(1,1); yt_frame(file_number,:) = yt + absolute_bbox(1,2);

%% Plotting the data

figure(1),imshow(i4), hold on, rectangle('Position', nextbboxref(fileCount,:),'EdgeColor','w'), plot(xh_frame(file_number,1),yh_frame(file_number,1),'r*'), plot(xt_frame(file_number,1),yt_frame(file_number,1),'g+'), hold off;
% figure(2), plot(yh_frame, xh_frame, 'r*'), axis([0 64 0 320]);

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
sprintf('Total Number of frames skipped : %d', skip_frame);