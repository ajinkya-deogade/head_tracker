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
j = 1;
if reply == 'Y'
    parse_number_of_frames = total_number_of_frames;
elseif reply == 'N'
    prompt2= 'Enter the number of frames to be extracted: ';
    parse_number_of_frames = input(prompt2);
end
%%  Loop to (1) extract & save frames, (2) mark co-ordinates of the point of interest (POI) and (3) save them in a file
coordinate = [];
figure(1), title('Please double click inside the box when done');
[I rect] = imcrop(read(vid_reader,1));
for i = 1:parse_number_of_frames
        vid_frame = read(vid_reader,i);
        vid_frame_crop = imcrop(vid_frame,rect);
        clc;
        disp(sprintf('Frame: %d',i));
        imwrite(vid_frame_crop, fullfile(pathname_jpg,[strcat(num2str(j,'%06d'),'_',name,'_', num2str(i)) ext_jpg]));
        j = j+1;
        close;
end
%% Release all the object handles and save data
close all;
clear all;
toc;