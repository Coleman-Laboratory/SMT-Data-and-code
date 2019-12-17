%%program to determine general rates of photobleaching based on overall
%%decay of fluorescence signal. User must first upload movie file, then set dt as the effect exposure
%%time of frames in movie.

clear all; close all;
[FileName, Pathname] = uigetfile('.tif'); %load image file
FileDTif=[Pathname, FileName];
prompt = 'What is the exposure time (seconds) of individual frames? ';
dt = input(prompt); %set exposure time


%% load stacked movie file --- DNA
h = waitbar(0,'Loading files...');

InfoImage=imfinfo(FileDTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage_DNA = zeros(nImage,mImage,NumberImages,'uint16');
 
TifLink = Tiff(FileDTif, 'r');
for i=1:NumberImages
   TifLink.setDirectory(i);
   temp = TifLink.read();
   maskedImage = temp; 
    FinalImage_DNA(:,:,i) = maskedImage;
    waitbar(i/NumberImages);
end
TifLink.close();
close(h);

%% calculate integrated intensity

first_ch = 6;
fig = 1;


imagesize = size(FinalImage_DNA, 3);
fluo = zeros(1,imagesize-first_ch+1);

for k = 1 : length(fluo)
    %fluo(k) = sum(sum(FinalImage_DNA(:,:,((k-1)*2+1)))); 
    fluo(k) = sum(sum(FinalImage_DNA(:,:,k+first_ch-1)));
end

%% fit to single exponential decay

t = dt * (first_ch:imagesize);

figure(fig);clf;
plot(t, fluo);

myfun = @(s) fluo - s(1) - s(2)*exp(-t/s(3));
f0 = [1e8, 4e8, 60];
options = optimset('Algorithm',{'levenberg-marquardt',0.000001});
s = lsqnonlin(myfun, f0, [], [], options);

yfit = s(1) + s(2)*exp(-t/s(3));

figure(fig); hold on;
plot(t, yfit, 'r')
hold off;
axis tight
xlabel('Time (sec)', 'FontSize', 16); 
ylabel('Fluo. Intensity', 'FontSize', 16);
set(gca, 'FontSize', 14);
title(['\tau_{bleach} = ', num2str(round(s(3))), ' sec'], 'FontSize', 18);



