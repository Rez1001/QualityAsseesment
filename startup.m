
clc
clear all
close all
% ---------------------- Read Files -------------------------

% PATH OF FOLDER
folder_orginal = 'D:\Retinal DataSets\final image my thesis\selected for image visualization\evaluation\Make random';
folder_enhance = 'D:\Retinal DataSets\final image my thesis\selected for image visualization\evaluation\CLAHE'; 

% YOU SHOULD CHOOSE TYPE OF YOUR ORGINAL AND ENHANCE  IMAGE
% BE CERFUL BOTH IMAGE SHOULD HAVE SAME SIZE ;)

filetype_orginal = '.jpg'; 
filetype_enhance = '.jpg'; 

files_orginal = dir(fullfile(folder_orginal, strcat('*',filetype_orginal))); 
files_enhance = dir(fullfile(folder_enhance, strcat('*',filetype_enhance))); 

 % number of files in the folder
num_files_orginal = numel(files_orginal);
num_files_enhance = numel(files_enhance);

images_orginal = cell(1,num_files_orginal); 
images_enhance = cell(1,num_files_enhance); 

for i = 1:num_files_orginal

    filename_orginal = fullfile(folder_orginal, files_orginal(i).name); 
    images_orginal{i} = imread(filename_orginal); 

end

for i = 1:num_files_enhance
     
    filename_enhance = fullfile(folder_enhance, files_enhance(i).name);  
    images_enhance{i} = imread(filename_enhance); 

end

% -------------------- USE INDEX YOU WANT ------------------

metric_names = {'Entropy Difference',"CII","Brisque Score",'Nique Score',"Brisque_score_orginal","Nique_score_orginal","psnr","ssim"};
results_cell=[metric_names];

if num_files_orginal == num_files_enhance

    disp('Numbr of Your orginal and enhance image is equal')

    for i = 1:num_files_orginal
    
         [Entropy_difference,CII,Brisque_score,Nique_score,Brisque_score_orginal,Nique_score_orginal,psnr_,ssim_] = MesurePerformance(images_orginal{i},images_enhance{i});
    
         metric_values = {Entropy_difference,CII,Brisque_score,Nique_score,Brisque_score_orginal,Nique_score_orginal,psnr_,ssim_};

         results_cell = [results_cell; metric_values];
         
    
    end
else
    disp('Numbr of Your orginal and enhance image is not equal')
end



% --------------------   FUNCTIONS  --------------------------- 

function [Entropy_difference,CII,Brisque_score,Nique_score,Brisque_score_orginal,Nique_score_orginal,ssim_,psnr_]=MesurePerformance(OriginalImage,OutputImage)

%Calculations of Entropy for both Images
Entropy_Original=entropy(OriginalImage);
Entropy_Output=entropy(OutputImage);
Entropy_difference=Entropy_Original-Entropy_Output;

%Contrast Improvement Index (CII)
CII=Calculate_CII(OriginalImage,OutputImage);

%Brisque score calculation
Brisque_score = brisque(OutputImage);

%Nique score calculation
Nique_score = niqe(OutputImage);

%Brisque score calculation
Brisque_score_orginal = brisque(OriginalImage);

%Nique score calculation
Nique_score_orginal = niqe(OriginalImage);

% psnr score calculation
psnr_=getPSNR(OriginalImage,OutputImage);

% ssim score calculation
ssim_=getMSSIM(OriginalImage,OutputImage);

end

% ------------- Contrast improvment index

function CII = Calculate_CII(Original_Image,Output_Image)

    kernel = ones(3)/9;
    
    % Compute sliding mean of original image
    slidingmean_Original = Original_Image;
    for dim = 1:ndims(Original_Image)
        slidingmean_Original = convn(slidingmean_Original, kernel, 'same');
    end
    
    % Compute sliding mean of output image
    slidingmean_Output_Image = Output_Image;
    for dim = 1:ndims(Output_Image)
        slidingmean_Output_Image = convn(slidingmean_Output_Image, kernel, 'same');
    end
    
    % Compute means and contrast improvement index (CII)
    mean_Original = mean(slidingmean_Original(:));
    mean_output_Image = mean(slidingmean_Output_Image(:));
    CII = mean_output_Image/mean_Original;
end

% ---------- SIIM ---------
function mssim = getMSSIM(frameReference,frameUnderTest)

C1 = 6.5025;
C2 = 58.5225;
frameReference=double(frameReference);
frameUnderTest=double(frameUnderTest);
frameReference_2=frameReference.^2;
frameUnderTest_2=frameUnderTest.^2;
frameReference_frameUnderTest=frameReference.*frameUnderTest;
mu1=imgaussfilt(frameReference,1.5);
mu2=imgaussfilt(frameUnderTest,1.5);
mu1_2=mu1.^2;
mu2_2=mu2.^2;
mu1_mu2=mu1.*mu2;
sigma1_2=imgaussfilt(frameReference_2,1.5);
sigma1_2=sigma1_2-mu1_2;
sigma2_2=imgaussfilt(frameUnderTest_2,1.5);
sigma2_2=sigma2_2-mu2_2;
sigma12=imgaussfilt(frameReference_frameUnderTest,1.5);
sigma12=sigma12-mu1_mu2;
t3 = ((2*mu1_mu2 + C1).*(2*sigma12 + C2));
t1 =((mu1_2 + mu2_2 + C1).*(sigma1_2 + sigma2_2 + C2));
ssim_map =  t3./t1;
mssim = mean2(ssim_map); mssim=mean(mssim(:));
end

% ------------------------ PSNR

function psnr=getPSNR(frameReference,frameUnderTest)

s1=double(frameReference-frameUnderTest).^2;
    
    s = sum(sum(s1)); 
    sse = s(:,:,1)+s(:,:,2)+s(:,:,3);
    if( sse <= 1e-10) 
        psnr=0;
    else
        mse  = sse / double(size(frameReference,1)*size(frameReference,2)*size(frameReference,3));
        psnr = 10.0 * log10((255 * 255) / mse);
    end
end

% -------------------------------



