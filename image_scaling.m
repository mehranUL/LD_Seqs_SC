N=256; % N<256 requires quantization  

sequences1 = transpose(net(sobolset(6),(N)));
sequences2 = VonDerCorput_6(N);
sequences3 = LFSR_Bulk(N)/N;
sequences4 = niederreiter2_generate(20, N, 2, 31);


%alpha values and image - FOREGROUND - F
[image, ~, alpha] = imread('your_alpha_image.png');
image = rgb2gray(image); %converting to grayscale
image = im2double(image); %converting image to double data format
alpha = im2double(alpha);

%BACKGROUND - B
background = imread('your_background_image.jpg');
background = rgb2gray(background);
background = im2double(background);

%Making the Image and Alpha matched to the Background in terms of dimensions
image = imresize(image, [size(background, 1) , size(background, 2)]);
alpha = imresize(alpha, [size(background, 1) , size(background, 2)]);

%--------------------------------------------------------------------------
%Conventional processing
reference_conventional = alpha.*image + (1-alpha).*background; % Image compositing formula
%--------------------------------------------------------------------------

%pre-allocation
SC_composite_image = zeros(size(background, 1), size(background, 2));

for i=1:1:size(background, 1)
    for j=1:1:size(background, 2)

        %Background in actual SC bitstreams
        %background_SC = random_number_BINOM(background(i,j),N);
        background_SC = number_source(background(i,j), N, sequences4(1, :));

        
        %In case of any overflow
        if image(i,j) > 1
            image(i,j) = 1;
        end
        %Image in actual SC bitstreams
        %image_SC = random_number_BINOM(abs(image(i,j)),N);
        image_SC = number_source(abs(image(i,j)), N, sequences4(1, :));

        
        %In case of any overflow
        if alpha(i,j) > 1
            alpha(i,j) = 1;
        end
        %Alpha in actual SC bitstreams
        %alpha_SC = random_number_BINOM(abs(alpha(i,j)),N);
        alpha_SC = number_source(abs(alpha(i,j)), N, sequences4(2, :));
         
        %Bit-by-bit processing using MUX
        temp = mux2_to_1_adder(background_SC, image_SC, alpha_SC);

        %Summing 1s in a bitstream
        SC_composite_image(i,j) = sum(temp);
        
    end
end

%PSNR value
p = psnr(uint8(SC_composite_image), uint8(N*reference_conventional));
s = ssim(uint8(SC_composite_image), uint8(N*reference_conventional));

imshow(uint8(SC_composite_image))
