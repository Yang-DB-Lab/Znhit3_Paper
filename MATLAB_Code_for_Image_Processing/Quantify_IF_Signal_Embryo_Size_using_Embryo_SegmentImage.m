%%%%%%%%%%%%%%%%%%%%%%%
%% Purpose of the function
%
%% The function takes 2 parameters:
%    * Protein IF image
%    * Binary picture for embryo segmentation
%% Make sure this fold contains 
%    * Prepared and checked binary image for embryo segmentation
%    * IF image a interested protein
% 
% This script uses embryo defined areas to calculate nuclear protein signals
% in these areas. The total number of pixel of one area as well as the 
% total IF signal in each embryo area will be calculated.
%
% A predefined segmentation picture of embryos (manually or by other
% IF channels) that deneate each embryo but without overlap is used to
% check what areas are included (with overlap) in one embryo.
% The total number of pixel of each area of an embryo as well as the 
% total IF signal in this area will be calculated.
%
% Finally, the table with 
%   * Total number of pixels in each embryonic area
%   * Total IF signal intensity of each area
%   * Embryo origin of each area
% will be saved as a .CSV file for future use.

% NOTE: MATLAB maximum name length is 63 characters; 
%       longer name will be truncated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Quantify_IF_signal_using_Embryo_SegmentImage(embryo_shape_filename, ...
    protein_IF_filename)

    %% Get the embryo "photo" for segmentation of each embryo
    %%% The embryo photo should be 8 bit binary (grayscale) image
    %%% with no holes
    embryo_photo = imread(embryo_shape_filename);
    if(max(max(embryo_photo) ) == 1) % binary photo
        embryo_bw = embryo_photo;
    else
        embryo_gray = im2gray(embryo_photo);
        ebmryo_bw_thresh = graythresh(embryo_gray)
        embryo_bw = imbinarize(embryo_gray, ebmryo_bw_thresh);
    end
    % erase small objects of small noise due to non-perfect segmentation
    embryo_bw = bwareaopen(embryo_bw, 50);
    % compare the binary mask picture to the original one
    imshowpair(embryo_photo, embryo_bw, "montage")
    % Label all embryo areas
    [embryo_area, embryo_num] = bwlabel(embryo_bw, 8);
    
    % show embryo labeling as a picture
    embryo_prop = regionprops(embryo_area, 'Centroid');
    imshow(embryo_bw)
    hold on
    for k = 1:numel(embryo_prop)
        c = embryo_prop(k).Centroid;
        text(c(1), c(2), sprintf('%d', k), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle');
    end
    hold off
    % Save it as a PDF
    saveas(gcf, embryo_shape_filename + "_Labeling.pdf")
    
    
    
    %% Get the IF image of the protein interested
    protein_IF = imread(protein_IF_filename);
    protein_gray = im2gray(protein_IF);
    imshowpair(protein_IF, protein_gray, "montage")
    
    
    
    %% The final calculation
    
    % The embryo_IF data contains 6 columns:
    %    * embryo_number
    %    * embryo_x_center
    %    * embryo_y_center
    %    * embryo_total_pixels_with_IF_signal
    %    * embryo_total_IF_intensity
    %    * embryo_mean_IF_intensity

    % Define the matrix to store the calculation results
    embryo_IF = zeros(embryo_num, 6);
    % first column is embryo area number
    embryo_IF(:,1) = 1:embryo_num;
    
    % Fill value for each embryo area
    for embryo_idx = 1:embryo_num
        c = embryo_prop(embryo_idx).Centroid;
        % 2nd column of embryo_area data is x_axis center of this area
        embryo_IF(embryo_idx, 2) = c(1);
        % 3rd column of embryo_IF data is x_axis center of this area
        embryo_IF(embryo_idx, 3) = c(2);
        % 4th column of embryo_IF data is total number of pixels in this area
        embryo_IF(embryo_idx, 4) = sum(sum(embryo_area == embryo_idx));
        % 5th column of embryo_area data is total nuclear protein IF inensity in this area
        area_mask = (embryo_area == embryo_idx);
        area_IF_data = protein_gray .* uint8(area_mask);
        embryo_IF(embryo_idx, 5) = sum(sum(area_IF_data));
    
        % get average IF signal intensity for this embryo
        embryo_IF(:, 6) = embryo_IF(:,5) ./ embryo_IF(:,4);
    end
    
    % Transfer the embryo_area matrix into table, then give it column names    
    embryo_IF_table = array2table(embryo_IF, ...
                                       "VariableNames", ...
                                       {'embryo_number', ...
                                       'embryo_x_center', ...
                                       'embryo_y_center', ... 
                                       'embryo_total_pixels_with_IF_signal', ...
                                       'embryo_total_IF_intensity',... 
                                       'embryo_mean_IF_intensity'});
    % write embryo summary into table
    writetable(embryo_IF_table, ...
              "Embryo_" + protein_IF_filename + "_IF_inensity.csv")

end
%% References
% 
% label the areas:
% https://blogs.mathworks.com/steve/2006/11/17/labeling-labeled-objects/
%
%
