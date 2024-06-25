%%%%%%%%%%%%%%%%%%%%%%%
%% Purpose of the function

%% The function takes 3 parameters:
%    * DAPI images
%    * Binary picture for embryo segmentation
%% Make sure this fold contains 
%    * the DAPI images 
%    * Prepared and checked binary image for embryo segmentation
% 
% This script uses DAPI to defined areas of nuclei.
% For each DAPI Z slice, the number of DAPI areas as well as their 
% positions will be noted in a matrix. 
%
% Another predefined segmentation picture of embryos (manually or by other
% IF channels) that deneate each embryo but without overlap is used to
% check what areas are included (with overlap) in one embryo
% The embryo original of each DAPI area will be aslo noted.
%
% Finally, Numbered different DAPI areas for each embryo will be counted
% This will be used to calculate how many blastomeres are in each embryo.
% This will be saved as a .CSV file for future use.
%
% The edges of embryos will be print out in a TIF image as well as edges 
% of different DAPI areas from all DAPI slices (The same DAPI area from
% adjecent DAPI image should indicate the nucleus of one same blastomere 
% and will only be printed for one time in the final edge image.
% 

% NOTE: MATLAB maximum name length is 63 characters; 
%       longer name will be truncated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Count_Nuclei_use_DAPI_and_Embryo_SegmentImage2(DAPI_filename_core, ...
    num_Z_slices, embryo_shape_filename)
    
    % get DAPI image size (height x width)
    DAPI_prime = imread(DAPI_filename_core + "1_20x.tif");
    DAPI_image_sz = size(DAPI_prime);
    % DAPI_image_sz
    % initilize a 3-d matrix to store all DAPI images
    DAPI_images = zeros(DAPI_image_sz(1), DAPI_image_sz(2), num_Z_slices);
    %% Get all DAPI images
    for i = 1:num_Z_slices

        DAPI = imread(DAPI_filename_core + i + "_20x.tif");
        DAPI_gray = im2gray(DAPI);
        bw_thresh = graythresh(DAPI_gray)
        DAPI_bw = imbinarize(DAPI_gray, max(bw_thresh, 0.1) );
    
        % Perform a morphological close operation on the image.
        % se = strel('disk',2);
        % DAPI_bw = imclose(DAPI_bw,se);
    
        %fill holes
        DAPI_bw = imfill(DAPI_bw,'holes');
    
        % compare the binary mask picture to the original one
        % imshowpair(DAPI_gray, DAPI_bw, "montage")
    
        % remove small areas which may be noise
        DAPI_bw_erased = bwareaopen(DAPI_bw, 30);

%         subplot(1,3,1), imshow(DAPI_gray)
%         subplot(1,3,2), imshow(DAPI_bw)
%         subplot(1,3,3), imshow(DAPI_bw_erased)
        % montage(DAPI_gray, DAPI_bw, DAPI_bw_erased)
        DAPI_images(:,:,i) = DAPI_bw_erased;
    end

    DAPI_projections = zeros(DAPI_image_sz(1), DAPI_image_sz(2), num_Z_slices);
    DAPI_projections_num_areas = zeros(1, num_Z_slices); % Note how many areas in each projection
    
    DAPI_num_areas = zeros(1, num_Z_slices);

    % Label all DAPI areas
    % use parameter 8 can reduce number of areas
    for i = 1:num_Z_slices
        DAPI_bw_erased = DAPI_images(:,:,i);
        [DAPI_images(:,:,i), DAPI_num_areas(i)] = bwlabel(DAPI_bw_erased, 8);
        % use DAPI_images directly to label DAPI areas
    end

    % For test only
%      for i=1:num_Z_slices
%          imwrite(DAPI_images(:,:,i), "BW_Frame1Channel5Slice" + i + "_20x.tif")
%      end

    % total_DAPI_projects = 1;

    DAPI_projections(:,:,1) = DAPI_images(:,:,1);
    DAPI_projections_num_areas(1) = DAPI_num_areas(1);
    

    % Check DAPI areas one by one and store image results to
    % "DAPI_projections"
    for i = 2:num_Z_slices
        % compare all DAPI areas of DAPI_image(:,:,i) 
        % to DAPI_image(:,:,i-1) and find new areas
        % Then update the new areas to DAPI_projections
        current_DAPI = DAPI_images(:,:,i);
        previous_DAPI = DAPI_images(:,:, i-1);
        

        for c=1:DAPI_num_areas(i) % iterate all DAPI areas of current DAPI bw

            current_DAPI_area_c = uint8(current_DAPI == c);
            current_DAPI_area_c_size = sum(sum(current_DAPI_area_c));
            % how many pixels in the DAPI area being examined
            % for current_DAPI bw image
                    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % compare with current Z projection images
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            handled_in_projection = false;% note weather one area is already handled in z projection
            zi=0;
            for z = 1:num_Z_slices

                DAPI_project_z = DAPI_projections(:,:,z); 
                if handled_in_projection == true
                    break; % already handled, do not need to check z-projections
                elseif max(max(DAPI_project_z)) == 0  & continue_to_next_z_projection == true% current projection is empty and the area not handled properly, add a new Z projection layer
                    DAPI_project_z(current_DAPI == c) = 1; % new layer built
                    % subplot(1,3,1), imshow(DAPI_projections(:,:,z))

                    % subplot(1,3,2), imshow(DAPI_project_z >0);
                    DAPI_projections(:,:,z) = DAPI_project_z;
                    DAPI_projections_num_areas(z) = DAPI_projections_num_areas(z) + 1;
                    handled_in_projection = true;
                    continue_to_next_z_projection = false;

                    subplot(1,2,1), imshow(DAPI_projections(:,:,1) >0)
                    subplot(1,2,2), imshow(DAPI_projections(:,:,2) >0)
                    break; % This area processed continue to next
                
                %%%%%%%%%%%%% The real job %%%%%%%%%
                else % check with current non-empty projection
                    [DAPI_projections(:,:,z), DAPI_projections_num_areas(z)] = bwlabel(DAPI_project_z>0, 8);
                    % label all areas in DAPI_project_z
                    continue_to_next_z_projection = false;
                    handled_in_projection = false;
                    for zi = 1:DAPI_projections_num_areas(z)
                        if continue_to_next_z_projection == true | handled_in_projection == true
                            break;
                        end

                        DAPI_projection_area_zi = uint8(DAPI_projections(:,:,z) == zi);
                        DAPI_projection_area_zi_size = sum(sum(DAPI_projection_area_zi));
                        % how many pixels in the DAPI area being examined
                        % for previous_DAPI bw image
                                
                        % How many overlapped pixels for these 2 DAPI areas
                        % being examined from "current_" and "previous_" DAPI
                        % bw images
                        num_overlap_zi = sum(sum(current_DAPI_area_c .* DAPI_projection_area_zi ));

                        %% A Z-projection DAPI area merges with the current checked DAPI area from current DAPI bw image
                        if num_overlap_zi >= 0.8* min(DAPI_projection_area_zi_size, current_DAPI_area_c_size) % Two nuclei merged together
                                                                               % check whehter it also overlaps previous DAPI image
                            num_overlap_pixel = 0; % current and previous DAPI area overlapped pixels
                            for p = 1:DAPI_num_areas(i-1) % all areas of previous DAPI bw
                                previous_DAPI_area_p = uint8(previous_DAPI == p);
                                previous_DAPI_area_p_size = sum(sum(previous_DAPI_area_p));
                                % how many pixels in the DAPI area being examined
                                % for previous_DAPI bw image
                
                                % How many overlapped pixels for these 2 DAPI areas
                                % being examined from "current_" and "previous_" DAPI
                                % bw images
                                num_overlap_pixel = sum(sum(current_DAPI_area_c .* previous_DAPI_area_p ));
                
                                if num_overlap_pixel >= 0.8 * min(previous_DAPI_area_p_size, current_DAPI_area_c_size) % and current DAPI area merges with previous DAPI
                                    if DAPI_projection_area_zi_size >= current_DAPI_area_c_size % Projection area bigger; do nothing
                                        continue_to_next_z_projection = false;
                                        handled_in_projection = true;
                                        break;
                                    else % Current DAPI projection is smaller, enlarge it
                                        DAPI_projection_z = DAPI_projections(:,:,z);       % same nucleus, enlarge its area in projection
                                        DAPI_projection_z(current_DAPI == c) = zi;
                                        DAPI_projections(:,:,z) = DAPI_projection_z;
                                        handled_in_projection = true;
                                        continue_to_next_z_projection = false;
                                        break;
                                    end
                                elseif num_overlap_pixel > 0 % Has small overlap with one previous DAPI area but not full overlap. New nucleus met, do sth on next z-projection
                                    continue_to_next_z_projection = true;
                                    handled_in_projection = false;
                                    break;
                                else % No overlap at all, continue to check next area in previous DAPI image
                                    continue;
                                end
                            end
                            if p== DAPI_num_areas(i-1) & num_overlap_pixel == 0 % completely new DAPI area compare to previous DAPI
                                continue_to_next_z_projection = true; % Although merge with current Z-projection ,it is a new nucleus.
                                handled_in_projection = false;
                                break;
                            end
                        elseif num_overlap_zi > 0 | continue_to_next_z_projection == true % overlap with current Z-projection but not merge. Go to next z-projection
                            continue_to_next_z_projection = true; % new nucleus met.
                            handled_in_projection = false;
                            break;
                        else
                            continue; % No overlap at all; continue to next Z-projection area
                        end
                    end
                    if zi == DAPI_projections_num_areas(z) & num_overlap_zi ==0 % continue_to_next_z_projection == false & handled_in_projection == false
                        % all zi checked but no overlap and not processed, means new area met
                        DAPI_projection_z = DAPI_projections(:,:,z); 
                        %%% IMPORTANT!!! Re-trieve it again, otherwise
                        %%% previous incorrect DAPI_projection_z for
                        %%% previous 'z' will be used !!!

                        DAPI_projection_z(current_DAPI == c) = zi+1; % add new area not exist before.
                        DAPI_projections(:,:,z) = DAPI_projection_z;
                        DAPI_projections_num_areas(z) = DAPI_projections_num_areas(z) + 1; % new area added; increase count.
                        handled_in_projection = true;
                        break; % This area handled. % Jump out of Z check.
                    end
                end
            end
        end
    end
    % for test
    for zz = 1:num_Z_slices
        imwrite(DAPI_projections(:,:,zz) > 0, "DAPI_Projection_" + zz + ".tif");
    end
    %%%%%%%%%%% Main work finished %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Begin label embryos and nuclei %%%%%%%%%%%%%%%%%%%%%%%
    % get embryo image
    embryo_image = imread(embryo_shape_filename);
    embryo_gray = im2gray(embryo_image);
    ebmryo_bw_thresh = graythresh(embryo_gray)
    embryo_bw = imbinarize(embryo_gray, ebmryo_bw_thresh);
    % erase small objects of small noise due to non-perfect segmentation
    embryo_bw = bwareaopen(embryo_bw, 50);

    % embryo edge
    embryo_edges = edge(embryo_bw); 
    
    %%%%%%%%%%%% Find out how many non-all-zero DAPI_projections
    total_non_zero_zz = 0;
    for zz = 1:num_Z_slices
        if max (max (DAPI_projections(:,:,zz))) > 0
            total_non_zero_zz = total_non_zero_zz + 1;
        end
    end

    % Define a series of nucleus colors
    text_colors = ["#FF0000","#00FFFF", "#D95319", "#FF00FF",... 
                   "#7E2F8E", "#EDB120", repelem("#000000" , 30) ];
    % "red",  "cyan","green", "magenta", 
    % finally just repeat "black" for enough times

    final_image = embryo_edges;

    %%%%%%%% Add different edges of DAPI projections on each layer
    %%%%%%%% to embryo_edges image.
    for zz = 1: total_non_zero_zz
        % get one DAPI projection
        one_DAPI_projection = DAPI_projections(:,:, zz) > 0;
        
        % find edges of DAPI projection area
        DAPI_projection_edges = edge(one_DAPI_projection);
        % "merge" DAPI edges into embryo edges image
        final_image = final_image + DAPI_projection_edges;
    end

    % imshowpair(1 - DAPI_projection_edges, 1 - embryo_edges, "montage")
    
    %% Technique to show all image together %%%%%%%%%%%%
    % "merge" DAPI_projection_edges and embryo_edges together
    % then reverse black white
    % imshow(1- (DAPI_projection_edges + embryo_edges) )
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    imshow(1 - final_image); % show embryo and DAPI area edges
    hold on    

    %% Label embryo numbers on picture
    [embryo_area, embryo_num] = bwlabel(embryo_bw, 8);
    
    % label embryo areas
    embryo_prop = regionprops(embryo_area, 'Centroid');
    for k = 1:embryo_num % numel(embryo_prop) also works
        c = embryo_prop(k).Centroid;
        embryo_k_bw = uint8(embryo_area == k);
        embryo_k_min_y = 0;
        for ky = 1: DAPI_image_sz(1) % embryo image has same size as DAPI
           if max(embryo_k_bw(ky,:))>0
               embryo_k_min_y = ky;
               break;
           end
        end
        text(c(1), embryo_k_min_y-5, sprintf('%d', k), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize',5, 'Color', 'b');
    end
    %% Label all DAPI areas
    for zz = 1:total_non_zero_zz
        [DAPI_area, DAPI_num_area] = bwlabel(DAPI_projections(:,:,zz), 8);
        % imshowpair(DAPI_projection, DAPI_projection_edges, "montage")
        % show embryo labeling as a picture
        DAPI_prop = regionprops(DAPI_area, 'Centroid');
        for k = 1:numel(DAPI_prop)
            c = DAPI_prop(k).Centroid;
            text(c(1), c(2), sprintf('%d', k), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 3, ... % use very small size
                'Color', text_colors(zz));
        end
   end

    hold off
    % Save it as a PDF
    saveas(gcf, embryo_shape_filename + "embryo_and_nucleus_Labeling.pdf")

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Quantification method %%%%%%%%%%%%%%%%%%
    
    %% The final calculation for summarization of each embryo
    % The embryo nucleus count summary contains 6 columns
    % total number of rows equals to toal number of nuclei 
    % from all DAPI projection layers

    % It contains the following 8 columns
    %    * DAPI area number
    %    * DAPI area layer number
    %    * DAPI area x_center value
    %    * DAPI area y_center value
    %    * DAPI area total number of pixels
    %    * DAPI area embryo orginal
    %    * x_center of the original embryo
    %    * y_center of the original embryo

    total_DAPI_areas = 0;
    for zz = 1:total_non_zero_zz
        zz_DAPI_area_num = max (max(DAPI_projections(:,:,zz)));
        total_DAPI_areas = total_DAPI_areas + zz_DAPI_area_num;
    end

    % Define the matrix to store the calculation results
    DAPI_area_data = zeros(total_DAPI_areas, 8);
    
    first_zz_row_idx = 0;
    for zz = 1:total_non_zero_zz
        [zz_DAPI_area, zz_DAPI_num] = bwlabel(DAPI_projections(:,:,zz), 8);
        % imshowpair(DAPI_projection, DAPI_projection_edges, "montage")
        % show embryo labeling as a picture
        DAPI_prop = regionprops(zz_DAPI_area, 'Centroid');
        
        % Fill value for each DAPI area
        for DAPI_idx = 1:zz_DAPI_num
            DAPI_area_data(first_zz_row_idx + DAPI_idx, 1) = DAPI_idx; % DAPI area number
            DAPI_area_data(first_zz_row_idx + DAPI_idx, 2) = zz; % DAPI area layer

            c = DAPI_prop(DAPI_idx).Centroid;
            % 3rd column of DAPI_area data is x_axis center of this area
            DAPI_area_data(first_zz_row_idx + DAPI_idx, 3) = c(1);
            % 4th column of DAPI_area data is x_axis center of this area
            DAPI_area_data(first_zz_row_idx + DAPI_idx, 4) = c(2);
            % 5th column of DAPI_area data is total number of pixels in this area
            DAPI_area_data(first_zz_row_idx + DAPI_idx, 5) = sum(sum(zz_DAPI_area == DAPI_idx));
        
            % get 6th column value (embryo original of DAPI area)
            for embryo_idx = 1:embryo_num
                embryo_idx_mask = (embryo_area == embryo_idx);
                area_mask = (zz_DAPI_area == DAPI_idx);
                if(sum(sum(uint8(area_mask) .* uint8(embryo_idx_mask) ) ) > 0) % DAPI_area has overlap with this embryo
                    DAPI_area_data(first_zz_row_idx + DAPI_idx, 6) = embryo_idx;   % This embryo contains this DAPI marked blastomere nuclear
                    % Put embryo center (x,y) value
                    embryo_center = embryo_prop(embryo_idx).Centroid;
                    DAPI_area_data(first_zz_row_idx + DAPI_idx, 7) = embryo_center(1);
                    DAPI_area_data(first_zz_row_idx + DAPI_idx, 8) = embryo_center(2);
                    break
                end
            end
        end
        first_zz_row_idx = first_zz_row_idx + zz_DAPI_num; % first row for next zz projection
    end
    
    % Transfer the DAPI_area_data matrix into table, then give it column names
    DAPI_area_table = array2table(DAPI_area_data, ...
                        'VariableNames',{'DAPI_area_number', ...
                                         'DAPI_area_layer', ...
                                         'x_center', ...
                                         'y_center', ...
                                         'total_pixels', ...
                                         'embryo_origin', ...
                                         'embryo_x_center', ...
                                         'embryo_y_center'});
    % Write the table to .csv file for future use
    writetable(DAPI_area_table, embryo_shape_filename + "_DAPI_area_and_embryo_origin.csv")
    % NOTE: 
    % If one DAPI area can not find its embryo origin (with embryo number 0),
    % that is because the DAPI area belongs to an embryo that is discarded
    % for analysis.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Another labeling method
    imshow(1 - final_image); % show embryo and DAPI area edges
    hold on
    %% Label embryo numbers on picture
    [embryo_area, embryo_num] = bwlabel(embryo_bw, 8);
    
    % label embryo areas
    embryo_prop = regionprops(embryo_area, 'Centroid');
    for embryo_idx = 1:embryo_num % numel(embryo_prop) also works
        c = embryo_prop(embryo_idx).Centroid;
        embryo_idx_bw = uint8(embryo_area == embryo_idx);
        embryo_idx_min_y = 0;
        for ky = 1: DAPI_image_sz(1) % embryo image has same size as DAPI
           if max(embryo_idx_bw(ky,:))>0
               embryo_idx_min_y = ky;
               break;
           end
        end
        text(c(1), embryo_idx_min_y-5, sprintf('%d', embryo_idx), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize',5, 'Color', 'b');

        % Begin label DAPI areas
        accumulated_DAPI_idx = 1;
        %% Label all DAPI areas based on its relationship with embryo
        for zz = 1:total_non_zero_zz
            [zz_DAPI_area, zz_DAPI_num_area] = bwlabel(DAPI_projections(:,:,zz), 8);
            % imshowpair(DAPI_projection, DAPI_projection_edges, "montage")
            % show embryo labeling as a picture
            DAPI_prop = regionprops(zz_DAPI_area, 'Centroid');
            for DAPI_idx = 1:numel(DAPI_prop)

                embryo_idx_mask = (embryo_area == embryo_idx);
                area_mask = (zz_DAPI_area == DAPI_idx);
                if(sum(sum(uint8(area_mask) .* uint8(embryo_idx_mask) ) ) > 0) % DAPI_area has overlap with this embryo
                    
                    DAPI_center = DAPI_prop(DAPI_idx).Centroid;
                    text(DAPI_center(1), DAPI_center(2), ...
                         sprintf('%d', accumulated_DAPI_idx), ...
                         'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'middle', ...
                         'FontSize', 3, ... % use very small size
                         'Color', text_colors(zz));
                    accumulated_DAPI_idx = accumulated_DAPI_idx + 1;
                end
            end
        end
    end
    hold off
    % Save it as a PDF
    saveas(gcf, embryo_shape_filename + "embryo_and_nucleus_Labeling2.pdf")


    