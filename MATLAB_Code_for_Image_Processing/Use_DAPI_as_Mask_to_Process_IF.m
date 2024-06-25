%%%%%%%%%%%%%%%%%%%%%%%
disp('Make sure this fold contains the DAPI file and OCT4 IF file with same base name')
disp('The files should be named like this: sample_DAPI.tif and sample_OCT4.tif')
prompt = "What is the sample name? ";
sample_name = input(prompt, "s")

%%%% USE DAPI as mask
DAPI = imread(sample_name + "_DAPI.tif");
DAPI_gray = im2gray(DAPI);
bw_thresh = graythresh(DAPI_gray)
DAPI_bw = imbinarize(DAPI_gray, bw_thresh/2);

% Perform a morphological close operation on the image.
se = strel('disk',2);
DAPI_bw = imclose(DAPI_bw,se);

%fill holes
DAPI_bw = imfill(DAPI_bw,'holes');

% compare the binary mask picture to the original one
imshowpair(DAPI_gray, DAPI_bw, "montage")

% remove polar body (which is smaller than nucleus) signal
DAPI_bw_erased = bwareaopen(DAPI_bw, 1000);
imshowpair(DAPI_gray, DAPI_bw_erased, "montage")

% get the original OCT IF picture
OCT4_IF = imread(sample_name + "_OCT4.tif");

% apply mask to remove non-sepcific signal
OCT4_IF_masked = bsxfun(@times, OCT4_IF, cast(DAPI_bw_erased, 'like', OCT4_IF));

% compare the original vs 'cleaned' OCT4 IF result
imshowpair(OCT4_IF, OCT4_IF_masked, "montage")

% save cleaned OCT4 IF result to a new .tif file. DONE!
imwrite(OCT4_IF_masked,sample_name + "_OCT4_masked.tif")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
