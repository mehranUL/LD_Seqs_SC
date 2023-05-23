% ICCAD
% UNDER REVIEW

clear all


N = 256;
sequences = transpose(net(sobolset(6),(N)));
%sequences = VonDerCorput_6(N);
%sequences = LFSR_Bulk(N)/N;
%sequences = niederreiter2_generate(20, N, 2, 31);

vidObj = VideoReader('kaplan.mp4');

%Creating a video
%REF: https://www.mathworks.com/matlabcentral/answers/280635-make-video-from-images
video = VideoWriter('new_sobol.mp4'); %create the video object
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
    alpha = double(alpha);
    alpha = alpha./255;

    background = double(rgb2gray(imread('zebra.jpg')));
    image = double(image);

    %matting
    %reference_conventional = alpha.*image + (1-alpha).*background; % Image compositing formula

    %pre-allocation
    SC_composite_image = zeros(size(background, 1), size(background, 2));

    for i=1:1:size(background, 1)
        for j=1:1:size(background, 2)

            %Background in actual SC bitstreams
            background_SC = number_source(background(i,j)./N, N, sequences(1, :));

            %Image in actual SC bitstreams
            image_SC = number_source(abs(image(i,j))./N, N, sequences(1, :));

            %Alpha in actual SC bitstreams
            alpha_SC = number_source(abs(alpha(i,j)), N, sequences(2, :));

            %Bit-by-bit processing using MUX
            temp = mux2_to_1_adder(background_SC, image_SC, alpha_SC);

            %Summing 1s in a bitstream
            SC_composite_image(i,j) = sum(temp);

        end
    end

    %PSNR value
    %p = psnr(uint8(SC_composite_image), uint8(N*reference_conventional));

    %figure
    %imshow(reference_conventional)
    writeVideo(video, uint8(SC_composite_image)); %write the image to file

end

close(video)