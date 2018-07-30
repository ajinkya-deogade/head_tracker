%% Input Files
close all;
clear all;
clc;
tic;
pathname_frames = uigetdir('./','Select the folder with frames');

%% Initialize Variables
dir_input = dir(fullfile(pathname_frames,'*.tiff'));
fileNames = {dir_input.name};
numFrames = numel(fileNames)
[path_folder, name_folder, ext_folder] = fileparts(pathname_frames);
marked_folder = fullfile(path_folder, strcat(name_folder,'_Cropped'));
mkdir(marked_folder);
cd(marked_folder);
j =1;
for i = 1:numFrames
    I = imread(fullfile(pathname_frames, fileNames{j}));
    I2 = imcrop(I,[165.5 52.5 65 208]);
    figure(1),
    imshow(I2);
    title(sprintf('Frame # %d',j));
    [path, name, ext] = fileparts(fileNames{j});
    imwrite(I2, fullfile(marked_folder,[strcat(name,'_crop') '.tiff']));
    j = j+1;
  end
close all;
toc;