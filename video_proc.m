% ICCAD
% UNDER REVIEW

clear all

vidObj = VideoReader('kaplan.mp4');

%Creating a video
%REF: https://www.mathworks.com/matlabcentral/answers/280635-make-video-from-images
video = VideoWriter('expected.mp4'); %create the video object
video.FrameRate = 30;
open(video); %open the file for writing


%vidObj.CurrentTime = 1;
% 
% while hasFrame(vidObj)
%     vidFrame = readFrame(vidObj);
%     imshow(vidFrame)
%     pause(1/vidObj.FrameRate);
% end

allFrames = read(vidObj);

[a, b, c, d] = size(allFrames);

for iter = 1:1:d

grayImage = allFrames(:,:,:,iter);

z = size(grayImage);
squeezed = reshape(grayImage, [540 960 3]);

grayImage = im2gray(squeezed);

%imshow(grayImage)

[row, column] = size(grayImage);

grayImageTemp = grayImage;


for i=1:1:row
    for j=1:1:column
        if      grayImage(i,j) == 129 || grayImage(i,j) == 130 || grayImage(i,j) == 131 || grayImage(i,j) == 132 || ...
                grayImage(i,j) == 133 || grayImage(i,j) == 134 || grayImage(i,j) == 135 || grayImage(i,j) == 136 || grayImage(i,j) == 137 || grayImage(i,j) == 138 || ...
                grayImage(i,j) == 139 || grayImage(i,j) == 140 ||  grayImage(i,j) == 141 || grayImage(i,j) == 142 ||  grayImage(i,j) == 143 || grayImage(i,j) == 144 || ...
                grayImage(i,j) == 145 || grayImage(i,j) == 146 || grayImage(i,j) == 147 || grayImage(i,j) == 148
            grayImageTemp(i,j) = 0;
        end
    end
end

%imshow(grayImageTemp)

image = grayImageTemp;

%alpha preparation
alpha = grayImageTemp;
for i=1:1:row
    for j=1:1:column
        if      grayImageTemp(i,j) ~= 0
            alpha(i,j) = 255;
        end
    end
end
%imshow(alpha)
alpha = alpha./256;


background = rgb2gray(imread('zebra.jpg'));

%matting
reference_conventional = alpha.*image + (1-alpha).*background; % Image compositing formula
%figure
%imshow(reference_conventional)
writeVideo(video, (reference_conventional)); %write the image to file

end

