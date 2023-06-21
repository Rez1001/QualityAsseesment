# QualityAsseesment
This is a simple code for QA between two series of images. only define your path and enjoy it. 

My code reads in pairs of images from two different folders, where each pair consists of an original image and an enhanced version of that image. It then computes several metrics to evaluate the quality of the enhancement process, including :

•	entropy difference
•	contrast improvement index (CII)
•	BRISQUE score
•	NIQE score
•	peak signal-to-noise ratio (PSNR)
•	structural similarity index (SSIM)

These metrics are calculated for each pair of images and stored in a matrix called 'results_cell'. 

