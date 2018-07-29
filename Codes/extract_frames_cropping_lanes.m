%% Input Files
close all;
clear all;
clc;
tic;
[filename_avi, pathname_avi, ext_avi] = uigetfile('*.avi','Select the video file to extract frames.');
pathname_jpg = uigetdir('./', 'Select the folder to save extracted image files.');
prompt1 = 'Do you want to extract all the frames (Y/N): ';
reply = input(prompt1,'s');
if isempty(reply)
    reply = 'Y';
end

%% Initialize Variables
ext_jpg='.jpg';
[path, name, ext_avi] = fileparts(filename_avi);
filename_jpg = fullfile(pathname_jpg,[name ext_jpg]);
vid_reader = VideoReader(strcat(pathname_avi, filename_avi));
lastFrame = read(vid_reader, inf);
total_number_of_frames = vid_reader.NumberOfFrames;

if reply == 'Y'
    parse_number_of_frames = total_number_of_frames;
elseif reply == 'N'
    prompt2= 'Enter the number of frames to be extracted: ';
    parse_number_of_frames = input(prompt2);
end
%%  Loop to (1) extract & save frames, (2) mark co-ordinates of the point of interest (POI) and (3) save them in a file
coordinate = [];
prompt2 = 'How many lanes do you want to extract: ';
I = read(vid_reader,1);
imshow(I);
num_lanes = input(prompt2);

for k = 1:(num_lanes)
    imshow(I);
    for l = 1:k-1
        rectangle('Position', rect_crop(l,:), 'EdgeColor','g');
    end
    messa = strcat('Draw box for lane ', num2str(k),' and double click when done.');
    annotation('textbox', [0.3, 0.8, .1, .1], 'String', messa);
    h = imrect(gca, [1, 1, 62, 162]);
    position = wait(h);
    rect_crop(k,:) = getPosition(h);
    close;
end
for  j = 1:num_lanes
    mkdir(pathname_jpg,strcat('Frame_', num2str(j)));
    for i = 1:parse_number_of_frames
        vid_frame = read(vid_reader,i);
        vid_frame_crop = imcrop(vid_frame, rect_crop(j,:));
        clc;
        disp(sprintf('Frame: %d',i));
        imwrite(vid_frame_crop, fullfile(pathname_jpg, strcat('Frame_', num2str(j)),[strcat(num2str(i,'%06d'),'_',name,'_', num2str(i)) ext_jpg]));
    end
end
%% Release all the object handles and save data
close all;
% clear all;
toc;