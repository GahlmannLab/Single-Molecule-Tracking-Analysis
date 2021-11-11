% Julian Rocha

%cellLocAnalysis

%Finds clusters in localization data and checks tracks position in relation
%to clusters

[dataFile, dataPath] = uigetfile({'*.mat';'*.*'},'Open all diffusion coefficients with tracks','MultiSelect', 'on');

[dataFile2, dataPath2] = uigetfile({'*.mat';'*.*'},'Open all fiberData files with cell outlines','MultiSelect', 'on');

tic

%align outlines to locs
totalLoc = struct;
totalDist = struct;
if iscell(dataFile2)
    for a = 1:length(dataFile2)
        load ([dataPath2,dataFile2{a}]);
        totalLoc(a).cell_info = fiberData;

    end
else
    load ([dataPath2,dataFile2]);
    totalLoc.cell_info = fiberData;
end

dist = cell(length(totalLoc),1);
rotLength = cell(length(totalLoc),1);
rotWidth = cell(length(totalLoc),1);
rotAxis = cell(length(totalLoc),1);
rotPts = cell(length(totalLoc),1);
rotTracks = cell(length(totalLoc),1);
totalDiffusion = [];
for s = 1:length(totalLoc)
    if iscell(dataFile)
            load ([dataPath,dataFile{s}]);
    else
        load ([dataPath,dataFile]);
    end
    
    if ~exist('tracks')
    tracks = tracks3;
    end
    
    totalDiffusion = vertcat(totalDiffusion,diffusionCoefficients);
    
    
    totInfo{s} = totalLoc(s).cell_info;
    cell_info_1dataSet = totalLoc(s).cell_info; % store one dataset information
    track_info_1dataSet = tracks; % store one dataset information JR
    data_length = length(cell_info_1dataSet);
    dist{s} = cell(data_length,1);
    dist_tracks = cell(data_length,1);
    ratio_dist = cell(data_length,1);
    ratio_dist_tracks = cell(data_length,1);
    outside_points = cell(data_length,1);
      

    
rotAxis{s} = cell(data_length,1);
rotPts{s} = cell(data_length,1);
rotTracks{s} = cell(data_length,1);
rotLength{s} = zeros(data_length,1);
rotWidth{s} = zeros(data_length,1);

    for m = 1:data_length % loop through one dataset cell by cell
      cell_info = cell_info_1dataSet(m);
        track_info = track_info_1dataSet(m);

        if mod(length(cell_info.meshPtsX),2) == 1
            cell_info.meshPtsX(end) = [];
            cell_info.meshPtsY(end) = [];
        end
        Num_points = length(cell_info.meshPtsX);
        cell_axis_mp = zeros(Num_points/2 - 1,3);
        cell_radius_vect_l = zeros(Num_points/2 - 1,2);
        cell_radius_vect_r = zeros(Num_points/2 - 1,2);
        cell_diameter_vect_s = zeros(Num_points/2 - 1,3);
        cell_diameter_vect_e = zeros(Num_points/2 - 1,3);
        cell_diameter_vect = zeros(Num_points/2 - 1,3);
        cell_axis_normal = zeros(Num_points/2 - 2,3);
        cell_diameter = zeros(Num_points/2 - 1,1);
        cell_radius = zeros(Num_points/2 - 1,1);
        theta=0:0.01:2*pi;
        
        cell_axis_mp_T = zeros(Num_points/2 - 1,3);
        cell_diameter_vect_s_T = zeros(Num_points/2 - 1,3);
        cell_diameter_vect_e_T = zeros(Num_points/2 - 1,3);
        cell_diameter_vect_T = zeros(Num_points/2 - 1,3);
        
        xLoc = cell_info.xLoc;
        yLoc = cell_info.yLoc;
        zLoc = cell_info.zLoc;

        xyz = [xLoc,yLoc,zLoc];
        xyz_proj = xyz;
        z_adjust = mean(zLoc); %%requires good filtering along z. AA
        Loc_axis_ind = zeros(size(xLoc)); % cell axis for each Loc
        Loc_diameter = zeros(size(xLoc)); % cell diameter for each Loc
        if length(xLoc) < 4
            continue
        else
            % start slice and end slice to select Locs for projection
            o = ceil(0.2*length(cell_axis_mp));  % start slice
            e = floor(0.8*length(cell_axis_mp));  % end slice
            
            % cell outlines from Oufti, add z coordinates
            cell_left_side_2d = [cell_info.meshPtsX(2:Num_points/2) cell_info.meshPtsY(2:Num_points/2)];
            cell_left_side_2d_new = zeros(size(cell_left_side_2d));
            cell_left_side_3d = [cell_left_side_2d z_adjust*ones(length(cell_left_side_2d(:,1)),1)];
            cell_right_side_2d = [cell_info.meshPtsX(Num_points/2 + 2:end) cell_info.meshPtsY(Num_points/2 + 2:end)];
            cell_right_side_2d = flipud(cell_right_side_2d);
            cell_right_side_2d_new = zeros(size(cell_right_side_2d));
            cell_right_side_3d = [cell_right_side_2d z_adjust*ones(length(cell_right_side_2d(:,1)),1)];
            
            %% Generate parameters for each slice before move the whole contour
            for i = 1:length(cell_left_side_3d(:,1))
                % Generate the center for each circle
                cell_axis_mp(i,:) = (cell_left_side_3d(i,:) + cell_right_side_3d(i,:))/2;
                % Generate the radius vector for each middle point
                cell_radius_vect_l(i,:) = cell_left_side_2d(i,:) - cell_axis_mp(i,1:2);
                cell_radius_vect_r(i,:) = cell_right_side_2d(i,:) - cell_axis_mp(i,1:2);
                % Starting point and end point of the diameter vector
                cell_diameter_vect_s(i,:) = cell_left_side_3d(i,:);
                cell_diameter_vect_e(i,:) = cell_right_side_3d(i,:);
                % Generate the vector of the diameter for each circle
                cell_diameter_vect(i,:) = cell_diameter_vect_s(i,:) - cell_diameter_vect_e(i,:);
                % Generate the length of diameter and radius
                cell_diameter(i) = sqrt((cell_left_side_3d(i,1) - cell_right_side_3d(i,1))^2 + (cell_left_side_3d(i,2) - cell_right_side_3d(i,2))^2 +...
                    (cell_left_side_3d(i,3) - cell_right_side_3d(i,3))^2);
                cell_radius(i) = 0.5*cell_diameter(i);
                % Extend the cell outline a little bit to generate corner coordinates for the polygon to check Loc in the following section
                cell_left_side_2d_new(i,:) = cell_left_side_2d(i,:) + 0.1*cell_radius_vect_l(i,:);
                cell_right_side_2d_new(i,:) = cell_right_side_2d(i,:) + 0.1*cell_radius_vect_r(i,:);
            end
            % Generate the normal vector for each circle. Indeed, this is not the normal vector as I thought, so I found a different way to
            % generate the normal vector.
            cell_axis_normal = diff(cell_axis_mp);
            points_cell = cell(length(cell_axis_normal),1);
            % Generate the normal vector of the slice that a Loc belongs to, based on 2D data
            % loop through points then slices
            %             figure; hold on; % plot polygons to check
            for p = 1:length(xLoc)
                xq = xLoc(p);
                yq = yLoc(p);
                for v = 1:length(cell_axis_normal) - 1
                    % create polygons to check Locs. Repeat first point so
                    % that all points in the shape are connected
                    xv = [cell_left_side_2d_new(v,1);cell_left_side_2d_new(v + 1,1);cell_right_side_2d_new(v + 1,1);cell_right_side_2d_new(v,1);cell_left_side_2d_new(v,1)];
                    yv = [cell_left_side_2d_new(v,2);cell_left_side_2d_new(v + 1,2);cell_right_side_2d_new(v + 1,2);cell_right_side_2d_new(v,2);cell_left_side_2d_new(v,2)];
                    in = inpolygon(xq,yq,xv,yv);
                    if in == 1
                        break
                    end
                end
                %                 plot(xv,yv);
                %                 plot(xq,yq,'ro');
                Loc_axis_ind(p) = v; %Writing which slice (v) the localization (p) resides in 
                Loc_diameter(p) = max(cell_diameter(v),cell_diameter(v + 1));
            end
            
            
 
            % Select Locs based on the start and end slices (subvolume the cell (use a loop for this section))
            for sb = e:-1:o+1
                volume_face_center_0 = cell_axis_mp(sb,:);  % start plane of the first subvolume
                volume_face_center_1 = cell_axis_mp(sb - 1,:);  % end plane of the first subvolume
                volume_Nvector_0 = cell_axis_normal(sb,:);
                volume_Nvector_1 = cell_axis_normal(sb - 1,:);
                volume_vector_0 = null(volume_Nvector_0); % Include two other unit vectors for its own orthogonal coordinate system
                volume_vector_1 = null(volume_Nvector_1); % Include two other unit vectors for its own orthogonal coordinate system
                
                % Generate corners (8) and center for the volume
                volumn_center = (volume_face_center_0 + volume_face_center_1)./2;
                
                volume_c1 = volume_face_center_0 + 500*volume_vector_0(:,1)' + 500*volume_vector_0(:,2)';
                volume_c2 = volume_face_center_0 + 500*volume_vector_0(:,1)' - 500*volume_vector_0(:,2)';
                volume_c3 = volume_face_center_0 - 500*volume_vector_0(:,1)' + 500*volume_vector_0(:,2)';
                volume_c4 = volume_face_center_0 - 500*volume_vector_0(:,1)' - 500*volume_vector_0(:,2)';
                volume_c5 = volume_face_center_1 + 500*volume_vector_1(:,1)' + 500*volume_vector_1(:,2)';
                volume_c6 = volume_face_center_1 + 500*volume_vector_1(:,1)' - 500*volume_vector_1(:,2)';
                volume_c7 = volume_face_center_1 - 500*volume_vector_1(:,1)' + 500*volume_vector_1(:,2)';
                volume_c8 = volume_face_center_1 - 500*volume_vector_1(:,1)' - 500*volume_vector_1(:,2)';
                
                corners_center = [volume_c1;volume_c2;volume_c3;volume_c4;volume_c5;volume_c6;volume_c7;volume_c8;volumn_center];
                
                x_tetra = corners_center(:,1);
                y_tetra = corners_center(:,2);
                z_tetra = corners_center(:,3);
                
                tri = delaunayTriangulation(x_tetra,y_tetra,z_tetra); % Generate delaunay triangulization
  
                tn = tsearchn(corners_center, tri, xyz_proj); % Determine which tetrahedron point is within
                IsInside = ~isnan(tn); % Convert to logical vector
                
                Loc_in = xyz_proj(IsInside,:);
                
                
                %% Project within Locs onto the start slice
                % normalize the normal vector of the start slice
                if ~isempty(Loc_in)
                    Nvect_o_slice = cell_axis_normal(sb - 1,:);
                    NNvect_o_slice = Nvect_o_slice/norm(Nvect_o_slice);
                    center_o_slice = cell_axis_mp(sb - 1,:); % center of the start slice
                    num_points_in = length(Loc_in(:,1)); % number of Locs between the start and end slices of the first subvolume
                    p_proj = zeros(num_points_in,3);
                    
                    % project within Locs to the start slice
                    for pp = 1:num_points_in
                        point_temp = Loc_in(pp,:);
                        p_proj(pp,:) = point_temp - dot(point_temp - center_o_slice, NNvect_o_slice) * NNvect_o_slice;
                    end
                    xyz_proj = [xyz_proj;p_proj];
                end
            end
            % Find the center (centroid) for projected Locs
            center_p_proj = [mean(p_proj(:,1)),mean(p_proj(:,2)),mean(p_proj(:,3))];
            
            %% Adjust cell outline according to the center of projected Locs and the center of the start slice
            % generate the translation vector
            T_vector = center_p_proj - center_o_slice;
            
            % translate cell outline
            cell_left_side_T = bsxfun(@plus,cell_left_side_3d,T_vector); % will this work?
            cell_right_side_T = bsxfun(@plus,cell_right_side_3d,T_vector); % will this work?
            
            % regenerate all parameters
            for i = 1:length(cell_left_side_T(:,1))
                % Generate the center for each circle
                cell_axis_mp_T(i,1) = (cell_left_side_T(i,1) + cell_right_side_T(i,1))/2;
                cell_axis_mp_T(i,2) = (cell_left_side_T(i,2) + cell_right_side_T(i,2))/2;
                cell_axis_mp_T(i,3) = (cell_left_side_T(i,3) + cell_right_side_T(i,3))/2;
                % Starting point and end point of the diameter vector
                cell_diameter_vect_s_T(i,:) = cell_left_side_T(i,:);
                cell_diameter_vect_e_T(i,:) = cell_right_side_T(i,:);
                % Generate the vector of the diameter for each circle
                cell_diameter_vect_T(i,:) = cell_diameter_vect_s_T(i,:) - cell_diameter_vect_e_T(i,:);
            end
            
            % Generate the normal vector for each circle
            cell_axis_normal_T = diff(cell_axis_mp_T);
            points_cell = cell(length(cell_axis_normal_T),1);
            
%                     % Plot all new circles (after translate)
%                     for j = 1:Num_points/2 - 2
%                         center = cell_axis_mp_T(j,:);
%                         radius = cell_radius(j);
%                         normal = cell_axis_normal_T(j,:);
%                         v = null(normal);
%                         points = repmat(center', 1, size(theta,2)) + radius*(v(:,1)*cos(theta) + v(:,2)*sin(theta));
%                         points_cell{j} = points;
%                         %     scatter3(center(1),center(2),center(3),'g');
%                         plot3(points(1,:), points(2,:), points(3,:),'g-');
%                     end




            midPt = ceil(length(cell_axis_mp_T)/2);
            Pt1 = midPt + ceil(length(cell_axis_mp_T)/4);
            Pt2 = midPt - ceil(length(cell_axis_mp_T)/4);

            if cell_axis_mp_T(Pt1,1) > cell_axis_mp_T(Pt2,1)
                theta = atan2((cell_axis_mp_T(Pt1,1)-cell_axis_mp_T(Pt2,1)),(cell_axis_mp_T(Pt1,2)-cell_axis_mp_T(Pt2,2)));
            else
                theta = atan2((cell_axis_mp_T(Pt2,1)-cell_axis_mp_T(Pt1,1)),(cell_axis_mp_T(Pt2,2)-cell_axis_mp_T(Pt1,2)));
            end
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                temp = cell_axis_mp_T(Pt2:Pt1,1:2);

                temp(:,1) = temp(:,1) - cell_axis_mp_T(midPt,1);
                temp(:,2) = temp(:,2) - cell_axis_mp_T(midPt,2);
                
                for h = 1:length(temp(:,1))
                  rotAxis{s}{m}(h,:) = (R*temp(h,:)')';
                end
                rotPts{s}{m} = xyz(:,:);

                rotPts{s}{m}(:,1) = rotPts{s}{m}(:,1) - cell_axis_mp_T(midPt,1);
                rotPts{s}{m}(:,2) = rotPts{s}{m}(:,2) - cell_axis_mp_T(midPt,2);
                rotPts{s}{m}(:,3) = rotPts{s}{m}(:,3) - cell_axis_mp_T(midPt,3);
                
                for h = 1:length(tracks{m})
                    rotTracks{s}{m}{h} = tracks{m}{h}(:,2:4);
                    rotTracks{s}{m}{h}(:,1) = rotTracks{s}{m}{h}(:,1) - cell_axis_mp_T(midPt,1);
                    rotTracks{s}{m}{h}(:,2) = rotTracks{s}{m}{h}(:,2) - cell_axis_mp_T(midPt,2);
                    rotTracks{s}{m}{h}(:,3) = rotTracks{s}{m}{h}(:,3) - cell_axis_mp_T(midPt,3);

                    for hh = 1:length(rotTracks{s}{m}{h})
                    rotTracks{s}{m}{h}(hh,1:2) = (R*rotTracks{s}{m}{h}(hh,1:2)')';
                    end
                end
                
                
                rotLength{s}(m) = abs(rotAxis{s}{m}(1,2)-rotAxis{s}{m}(end,2));
                rotWidth{s}(m) = cell_radius(ceil(length(temp(:,1))/2));
                for hh = 1:length(rotPts{s}{m}(:,1))
                    rotPts{s}{m}(hh,1:2) = (R*rotPts{s}{m}(hh,1:2)')';
                   
                end
                temp2 = rotPts{s}{m};
                rotPts{s}{m} = temp2(temp2(:,2) < max(rotAxis{s}{m}(:,2)) & temp2(:,2) > min(rotAxis{s}{m}(:,2)),:);
            
            % Calculate the distance between within Locs to the translated cell axis vectors
            num_points_1cell = length(xLoc);
            num_slice = length(cell_diameter_vect_T);
            
            for j = 1:num_points_1cell
                vector_sb = cell_axis_normal_T(Loc_axis_ind(j),:);
                % Calculate the distance of this Loc to the cell axis of the subvolume it belongs to
                b2 = xyz(j,:) - cell_axis_mp_T(Loc_axis_ind(j),:);
                dist_temp = norm(cross(vector_sb,b2)) / norm(vector_sb);
                point_temp = xyz(j,:);
                dist{s}{m}(j) = dist_temp;
                
                
                % Calculate the ratio of the distance to the cell radius at that Loc
                
                ratio_dist{m}(j) = dist_temp/Loc_diameter(j);
%                 if ratio_dist{m}(j) <= 1.1 % check if the distance is smaller than the radius,if not, print warning.
                if ratio_dist{m}(j) <= 1.25 % check if the distance is smaller than the radius,if not, print warning.
                else
                    outside_points{m}(j,:) = point_temp;
                    out_inf = ['FOV ', num2str(s), ' Cell ', num2str(m), ' Loc ', num2str(j), ' is outside the subvolume, Ratio: ', num2str(ratio_dist{m}(j))];
                    display(out_inf);
                    ratio_dist{m}(j) = 1;
                    %                 xout = outside_points{m}(j,1);
                    %                 yout = outside_points{m}(j,2);
                    %                 zout = outside_points{m}(j,3);
                    %                 plot3(xout,yout,zout,'go');
                end
            end
            
            dist_temp_2 = dist{s}{m}';
            dist_temp_3 = ratio_dist{m}';
            if isempty(track_info{:})
            else
                num_tracks_1cell = length(track_info{:});
                for c = 1:num_tracks_1cell % loop through tracks in cell m
                    num_Loc_1track = length(track_info{1}{c}(:,1));
                    for d = 1:num_Loc_1track % loop through Loc in each track c in cell m
                        xyz_track = track_info{1}{c}(d,2:4);
                        if isempty(strmatch(xyz_track,xyz))
                            track_info{1}{c}(d,6) = nan;
                            track_info{1}{c}(d,7) = nan;
                        else
                            track_info{1}{c}(d,6) = dist_temp_2(strmatch(xyz_track,xyz));
                            track_info{1}{c}(d,7) = dist_temp_3(strmatch(xyz_track,xyz));
                        end
                    end                    
                end
                totalTrackDD(s).tracks{m} = track_info{:};
            end
        end
    end
    totalDist(s).dist = dist{s};
    totalDist(s).dist_ratio = ratio_dist;
    totalDist(s).outside_points = outside_points;
    
    clear tracks tracks3
end
% save('Tracks_DC_distance_ratio.mat','totalDist','totalTrackDD'); % totalTrackDD contains tracks, diffusion coefficients, distance and distance ratio

%Calculate the scaling factors for each individual cell. Scale by the
%average lengths and widths for all FOV
rotScaleY = [];
rotScaleX = [];
meanScaleY = mean(cellfun(@mean, rotLength));
meanScaleX = mean(cellfun(@mean, rotWidth));
for s = 1:length(rotPts)
%         maxScaleY = mean(rotLength{s});
%         maxScaleX = mean(rotWidth{s});

    rotScaleY{s} = (meanScaleY.*ones(length(rotLength{s}),1))./rotLength{s};
    rotScaleX{s} = (meanScaleX.*ones(length(rotWidth{s}),1))./rotWidth{s};
end

%Apply the scaling factor in each dimension for each cell
rotPtsScale = [];
rotTracksScale = [];
for s = 1:length(rotPts)
    for m = 1:length(rotPts{s})
        if ~isempty(rotPts{s}{m})
        rotPtsScale{s}{m}(:,2) = rotPts{s}{m}(:,2)*rotScaleY{s}(m);
        rotPtsScale{s}{m}(:,1) = rotPts{s}{m}(:,1)*rotScaleX{s}(m);
        rotPtsScale{s}{m}(:,3) = rotPts{s}{m}(:,3)*rotScaleX{s}(m);
        end
        if ~isempty(rotTracks{s}{m})
            for h = 1:length(rotTracks{s}{m})
                rotTracksScale{s}{m}{h}(:,2) = rotTracks{s}{m}{h}(:,2)*rotScaleY{s}(m);
                rotTracksScale{s}{m}{h}(:,1) = rotTracks{s}{m}{h}(:,1)*rotScaleX{s}(m);
                rotTracksScale{s}{m}{h}(:,3) = rotTracks{s}{m}{h}(:,3)*rotScaleX{s}(m);
            end
        else
            rotTracksScale{s}{m} = [];
        end
    end
end

%Plot top view of 50% of the cell overlay. Cells are aligned to be
%vertical.
figure;
hold on;
%totRotTracksScale and totRotPts Scale is a cell array with information
%from all cells in all FOVs
totRotTracksScale = [];
totRotPtsScale = [];
totRotAxis = [];
totRotIdx = [];
idxFOV = [];
idxCell = [];
fThresh = 150;
for s = 1:length(rotPts)
%     figure; hold on;
    for m = 1:length(rotPts{s})
            if ~isempty(rotAxis{s}{m}) && rotAxis{s}{m}(1,1) > -fThresh && rotAxis{s}{m}(end,1) < fThresh
%                             plot(rotPtsScale{s}{m}(:,1),rotPtsScale{s}{m}(:,2),'.');
                            plot(rotPtsScale{s}{m}(:,1),rotPtsScale{s}{m}(:,3),'.');
                totRotAxis = vertcat(totRotAxis,rotAxis{s}{m});
                totRotPtsScale = vertcat(totRotPtsScale,rotPtsScale{s}{m});
                totRotIdx = vertcat(totRotIdx,ones(length(rotTracksScale{s}{m}'),1));
                totRotTracksScale = vertcat(totRotTracksScale,rotTracksScale{s}{m}');
                idxFOV = [idxFOV;s*ones(length(rotTracksScale{s}{m}'),1)];
                idxCell = [idxCell;m*ones(length(rotTracksScale{s}{m}'),1)];
                
            else
%                 totRotIdx = vertcat(totRotIdx,zeros(length(rotTracksScale{s}{m}'),1));
%                 idxFOV = [idxFOV;s*ones(length(rotTracksScale{s}{m}'),1)];
%                 idxCell = [idxCell;m*ones(length(rotTracksScale{s}{m}'),1)];

            end
%         title([num2str(s) '_' num2str(m)]);
%         axis equal;
%         pause
    end
%      title([num2str(s) '_' num2str(m)]);
%         axis equal;
%         pause
end
axis equal;

%set lim to crop out localizations far away from the central axis
Lim = 1000;
cntrs{1} = -1000:20:1000;
cntrs{2} = -1000:20:1000;
xticklabels = cntrs{1}(1):200:cntrs{1}(end);
xticks = 0:10:100;
yticklabels = cntrs{2}(1):200:cntrs{2}(end);
yticks = 0:10:100;


% %Plot side view of 50% of the cell overlay, 1 cell at a time. Also plot
% %mean of localizations. Used to verify good overlay.
% figure;
% for s = 1:length(rotPts)
%     for m = 1:length(rotPts{s})
%         if ~isempty(rotAxis{s}{m}) && rotAxis{s}{m}(1,1) > -fThresh && rotAxis{s}{m}(end,1) < fThresh
% %                 hold on;
%                 plot(rotPtsScale{s}{m}(:,1),rotPtsScale{s}{m}(:,3),'.b');
% %                 plot(mean(rotPtsScale{s}{m}(:,1)),mean(rotPtsScale{s}{m}(:,3)),'.r','MarkerSize',30);
% %                 hold off
%         end
%         axis equal;
%         
%         pause;
%     end
% end
% xlim([-Lim Lim]);
% ylim([-Lim Lim]);


%Plot slow vs fast tracks side view of 50% of the cell overlay
slowThresh = 0.15;
figure;
hold on;
totRotTrackPts = [];
totRotTracksSlow =[];
totRotTracksFast =[];
totalDiffusionCrop = totalDiffusion(logical(totRotIdx));
tempSlowIdx = [];
for s = 1:length(totRotTracksScale)
    if totalDiffusionCrop(s) < slowThresh & totRotIdx(s) == 1
        plot(totRotTracksScale{s}(:,1),totRotTracksScale{s}(:,3),'.b');
        totRotTracksSlow = vertcat(totRotTracksSlow,totRotTracksScale{s});
        tempSlowIdx = [tempSlowIdx;s];
    else
        plot(totRotTracksScale{s}(:,1),totRotTracksScale{s}(:,3),'.r');
        totRotTracksFast = vertcat(totRotTracksFast,totRotTracksScale{s});
    end
%             plot(totRotTracksScale{s}(:,1),totRotTracksScale{s}(:,3),'.');
            totRotTrackPts = vertcat(totRotTrackPts,totRotTracksScale{s});
end
axis equal;

%For troubleshooting slow trajectory clusters
slowTrackCell = idxCell(tempSlowIdx);
tempIdx = 1;
figure;
for s = 1:numCells
    while slowTrackCell(tempIdx) == s
        plot(totRotTracksScale{tempSlowIdx(tempIdx)}(:,1),totRotTracksScale{tempSlowIdx(tempIdx)}(:,3),'.b');
        hold on;
        axis equal;
        title(num2str(s));
%         pause;
        tempIdx = tempIdx + 1;
        if tempIdx > length(slowTrackCell)
            break;
        end
    end
    if s > max(slowTrackCell)
        break;
    end
%     hold off;
end
close;

%Crop out localizations in trajectories far away from the central axis
totRotPtsScale = totRotPtsScale(totRotPtsScale(:,1)>-Lim & totRotPtsScale(:,1)<Lim & totRotPtsScale(:,3)>-Lim & totRotPtsScale(:,3)<Lim,:);
totRotTrackPts = totRotTrackPts(totRotTrackPts(:,1)>-Lim & totRotTrackPts(:,1)<Lim & totRotTrackPts(:,3)>-Lim & totRotTrackPts(:,3)<Lim,:);
totRotTracksSlow = totRotTracksSlow(totRotTracksSlow(:,1)>-Lim & totRotTracksSlow(:,1)<Lim & totRotTracksSlow(:,3)>-Lim & totRotTracksSlow(:,3)<Lim,:);
totRotTracksFast = totRotTracksFast(totRotTracksFast(:,1)>-Lim & totRotTracksFast(:,1)<Lim & totRotTracksFast(:,3)>-Lim & totRotTracksFast(:,3)<Lim,:);

%plot side view of 50% of the cell overlay, 3D histogram
% figure; temp3 = hist3([totRotPtsScale(:,3) totRotPtsScale(:,1)],[100 100]); axis equal; pcolor(temp3); axis equal;
figure; tempSideLoc = hist3([totRotPtsScale(:,3) totRotPtsScale(:,1)],cntrs); ...
    imagesc(flipud(tempSideLoc)); set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);...
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);axis equal; ...
    colormap gray; title('All Locs'); xlabel('nm'); ylabel('nm');


%Plot tracks side view of 50% of the cell overlay, 3D histogram
figure; temp3 = hist3([totRotTrackPts(:,3) totRotTrackPts(:,1)],[100 100]); imagesc(flipud(temp3)); ...
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);...
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);axis equal; ...
    colormap gray; title('All Track Locs'); xlabel('nm'); ylabel('nm');
figure; temp4 = hist3([totRotTracksSlow(:,3) totRotTracksSlow(:,1)],[100 100]); imagesc(flipud(temp4));...
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);...
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);axis equal; ...
    colormap gray; title('Slow Track Locs'); xlabel('nm'); ylabel('nm');
figure; temp5 = hist3([totRotTracksFast(:,3) totRotTracksFast(:,1)],[100 100]); imagesc(flipud(temp5));...
    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);...
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels);axis equal; ...
    colormap gray; title('Fast Track Locs'); xlabel('nm'); ylabel('nm');

% %Plot 2D histograms of the side views of the tracks in previous figures.
% temp3_row = sum(temp3,2); temp3_col = sum(temp3,1); 
% figure; plot(1:length(temp3),temp3_row); figure; plot(1:length(temp3),temp3_col);
% temp4_row = sum(temp4,2); temp4_col = sum(temp4,1); 
% figure; plot(1:length(temp4),temp4_row); figure; plot(1:length(temp4),temp4_col);
% temp5_row = sum(temp5,2); temp5_col = sum(temp5,1); 
% figure; plot(1:length(temp5),temp5_row); figure; plot(1:length(temp5),temp5_col);

%Calculate the scaled (to the mean width) distance of the track points to the central axis
totRotTrackScaleDist = sqrt(totRotTrackPts(:,1).^2 + totRotTrackPts(:,3).^2);
totRotTrackSlowScaleDist = sqrt(totRotTracksSlow(:,1).^2 + totRotTracksSlow(:,3).^2);
totRotTrackFastScaleDist = sqrt(totRotTracksFast(:,1).^2 + totRotTracksFast(:,3).^2);

%Plot the count histograms scaled (to the mean width) distance of the track points to the central axis
edges = 0:10:Lim;
binspacing = edges(2);
figure; hold on;
title('Distance of Track Locs to Central Axis');
histogram_counts = histcounts(totRotTrackFastScaleDist,edges);
pdf_distribution = histogram_counts;
bar1 = bar((edges(2:end)-binspacing/2),pdf_distribution,'r');
histogram_counts = histcounts(totRotTrackSlowScaleDist,edges);
pdf_distribution = histogram_counts;
bar1 = bar((edges(2:end)-binspacing/2),pdf_distribution,'b');
line([meanScaleX, meanScaleX], ylim, 'LineWidth', 2, 'Color', 'k','LineStyle','--');
legend('Fast Tracks','Slow Tracks (<0.15)','Mean Cell Radius');
xlabel('Distance (nm)');
ylabel('Counts');

%Plot the normalized PDF histograms scaled (to the mean width) distance of the track points to the central axis
edges = 0:10:Lim;
binspacing = edges(2);
figure; hold on;
title('PDF Normalized Distance of Track Locs to Central Axis');
histogram_counts = histcounts(totRotTrackFastScaleDist,edges);
pdf_distribution = histogram_counts./(sum(histogram_counts));
bar1 = bar((edges(2:end)-binspacing/2),pdf_distribution,'r');
histogram_counts = histcounts(totRotTrackSlowScaleDist,edges);
pdf_distribution = histogram_counts./(sum(histogram_counts));
bar1 = bar((edges(2:end)-binspacing/2),pdf_distribution,'b');
line([meanScaleX, meanScaleX], ylim, 'LineWidth', 2, 'Color', 'k','LineStyle','--');
legend('Fast Tracks','Slow Tracks (<0.15)', 'Mean Cell Radius');
xlabel('Distance (nm)');
ylabel('Probability');


% %plot cell axis overlay
% figure;
% hold on;
% for s = 1:length(rotAxis)
% for m = 1:length(rotAxis{s})
% plot(rotAxis{s}{m}(:,1),rotAxis{s}{m}(:,2));
% axis equal;
% pause
% end
% end
% axis equal;


%% Old code, cluster localization analysis


totalDistance = [];
totalRatio = [];
for j = 1:length(totalDist)
    for k = 1:length(totalDist(j).dist)
% for j = 1:length(totalDist.dist)
%     distance = cell2mat(totalDist.dist(j))';
    distance = (totalDist(j).dist{k})';
    totalDistance = vertcat(totalDistance,distance);
    ratio = (totalDist(j).dist_ratio{k})';
    totalRatio = vertcat(totalRatio,ratio);
    end
end

% dbstop if error

time_per_frame = 25;


% diffThresh = 0.25;
diffThresh = 0.15;
% distThresh = 100;
distThresh = 150;


totalDiffusion = [];
totalTrackLength = [];
numCellsTotal = 0;
tracksTotal = [];
cellIdx = [];
idx = 0;
cellOutlines = [];
cellLocs = [];
if iscell(dataFile)
    for a = 1:length(dataFile)
        load ([dataPath,dataFile{a}]);
        load ([dataPath2,dataFile2{a}]);
        totalDiffusion = vertcat(totalDiffusion,diffusionCoefficients);
        totalTrackLength = vertcat(totalTrackLength,totTrackLength);
        numCellsTotal = numCellsTotal + numCells;
        
        for b = 1:length(fiberData)
        idx = idx + 1;
        cellIdx = vertcat(cellIdx,idx*ones(length(tracks{b}),1));
        tracksTotal = vertcat(tracksTotal,tracks{b});
        tempMesh = [];
        tempMesh{1} = [fiberData(b).meshPtsX, fiberData(b).meshPtsY];
        cellOutlines = vertcat(cellOutlines, tempMesh);
        tempLocs = [];
        tempLocs{1} = [fiberData(b).xLoc, fiberData(b).yLoc, fiberData(b).zLoc];
        cellLocs = vertcat(cellLocs,tempLocs);
        end
    end
else
    load ([dataPath,dataFile]);
    totalDiffusion = vertcat(totalDiffusion,diffusionCoefficients);
    totalTrackLength = totTrackLength;
    numCellsTotal = numCellsTotal + numCells;
    for b = 1:length(fiberData)
        idx = idx + 1;
        tracksTotal = vertcat(tracksTotal,tracks{b});
        cellIdx = idx*ones(length(tracks{b}),1);
        cellLocs{b} = [fiberData(b).xLoc, fiberData(b).yLoc, fiberData(b).zLoc];
        cellOutlines{b} = [fiberData(b).meshPtsX, fiberData(b).meshPtsY];
    end
end

clear totTrackLength


slowTrackIDX = zeros(length(totalDiffusion),1);
clustIDX = cell(numCellsTotal,1);
clustCOM = cell(numCellsTotal,1);
clustPtsIDX = cell(numCellsTotal,1);
for c = 1:numCellsTotal    
    tempTracks = tracksTotal(cellIdx == c);
    tempDiff = totalDiffusion(cellIdx == c);
    
% % %     epsilon = 50;
% % %     MinPts = 15;
% %     epsilon = 75;
% %     MinPts = 20;
%     epsilon = 50;
%     MinPts = 10;
    epsilon = 75;
    MinPts = 15;

    X = [cellLocs{c}(:,1),cellLocs{c}(:,2),cellLocs{c}(:,3)];
    if length(X(:,1)) == 1
        continue;
    end
    
    [IDX, isnoise]=DBSCAN(X,epsilon,MinPts);
    IDX2 = logical(IDX);
    
    if c == 1
            clustIDX{c} = IDX ;
            maxIDX = 0;
    else
            if sum(clustIDX{c-1}) > 0
            maxIDX = max(clustIDX{c-1});
            end
            clustIDX{c} = (IDX~=0).*(IDX + maxIDX);

    end
%     clustIDX{c} = IDX2;
    
    tempIDX = [];
    clustCOM{c} = zeros(max(IDX),3);
    for j = 1:max(IDX)
        tempIDX = logical(IDX == j);
        clustCOM{c}(j,:) = mean(cellLocs{c}(tempIDX,:)); 
    end
    
%     figure;
%     plot(cellOutlines{c}(:,1),cellOutlines{c}(:,2),'-g');
%     hold on;
%     plot3(cellLocs{c}(:,1),cellLocs{c}(:,2),cellLocs{c}(:,3),'.r');
%     axis equal;
% 
%     figure;
%     plot(cellOutlines{c}(:,1),cellOutlines{c}(:,2),'-g');
%     hold on;
%     plot3(cellLocs{c}(IDX2,1),cellLocs{c}(IDX2,2),cellLocs{c}(IDX2,3),'.b');
%     axis equal
    
    idxTrackTemp = zeros(length(tempTracks),1);
    idxTrackTemp(tempDiff <= diffThresh) = 1;
    idxTrackTempNum = zeros(length(tempTracks),1);
    if ~isempty(clustCOM{c})
        for e = 1:length(clustCOM{c}(:,1))
            for f = 1:length(tempTracks)
                for g = 1:length(tempTracks{f})
                    distance = pdist2(tempTracks{f}(g,2:4),clustCOM{c}(e,:));
                    if distance <= distThresh 
%                            idxTrackTemp(f,1) = 2;
%                            clustPtsIDX{c}{f}(g,1) = 1;
                            if idxTrackTemp(f,1) == 1
                            idxTrackTemp(f,1) = 2;
                            end
%                             clustPtsIDX{c}{f}(g,1) = 1;
                    end
                end
            end
        end
        
        %slowTrackIDX = 1 means track diffusion coefficient is below
        %threshold
        %slowTrackIDX = 2 means track diffusion coefficient is below
        %threshold and is close to a cluster of localizations
        slowTrackIDX(cellIdx == c) = idxTrackTemp;
        
    end
     
    z = 1;
    
    c
end
toc;

%slowTrackClust the percentage of all slow trajectories close to a cluster
slowTrackClust = sum(slowTrackIDX == 2)/sum(slowTrackIDX > 0).*100


%slowClustPerc is the percentage of all clustered points that are part of a
%slow trajectory
clustIDXtot = [];
for j = 1:numCellsTotal
    clustIDXtot = vertcat(clustIDXtot,clustIDX{j});
end
slowClustPerc = sum((slowTrackIDX == 2).*totalTrackLength)/sum(clustIDXtot ~= 0).*100


% figure; 
% edges = 0:0.05:12;
% binSpacing = edges(2)-edges(1);
% totalDiffusionCounts = histcounts(totalDiffusion,(edges+binSpacing/2));
% totalDiffusionCounts = totalDiffusionCounts./sum(totalDiffusionCounts);
% bar(edges(2:end),totalDiffusionCounts);
% % hold on; histogram(totalDiffusion(trackIDX == 1),(edges+binSpacing/2)); 
% % hold on; histogram(totalDiffusion(trackIDX == 2),(edges+binSpacing/2));
% totalBindCounts = histcounts(totalDiffusion(trackIDXBind > 1),(edges+binSpacing/2));
% totalBindCounts = totalBindCounts./sum(totalBindCounts);
% hold on; bar(edges(2:end),totalBindCounts,'r');

%plot distance to cell axis for non-clustered points
total_num_dist = length(totalDistance);
figure;
hold on;
title(['Distance to the cell axis']);
xlabel('Distance (nm)');
ylabel('Probability Density');
edges = 0:10:max(totalDistance);
binSpacing = edges(2)-edges(1);
histogram_counts = histcounts(totalDistance(clustIDXtot == 0),(edges+binSpacing/2));
pdf_distribution = histogram_counts;
% pdf_distribution = histogram_counts./(sum(histogram_counts));
bar1 = bar(edges(2:end),pdf_distribution,'r');

%plot distance to cell axis for cluster points
% totalDistanceClust = totalDistance(clustIDXtot == 1);
totalDistanceClust = totalDistance(clustIDXtot ~= 0);
total_num_dist_clust = length(totalDistanceClust);
% figure;
hold on;
% title(['Distance to the cell axis']);
% xlabel('Distance (nm)');
% ylabel('Probability Density');
% edges = 0:10:max(totalDistanceClust);
% binSpacing = edges(2)-edges(1);
histogram_counts = histcounts(totalDistanceClust,(edges+binSpacing/2));
pdf_distribution = histogram_counts;
% pdf_distribution = histogram_counts./(sum(histogram_counts));
bar1 = bar(edges(2:end),pdf_distribution,'b');
xlim([0 Lim])

%plot stationary vs fast tracks 10/25
totalTracks = [];
e = 1;
for a = 1:length(totalTrackDD)
for b = 1:length(totalTrackDD(a).tracks)
    for c = 1:length(totalTrackDD(a).tracks{b})        
        totalTracks{e,1} = totalTrackDD(a).tracks{b}{c};
        e = e+1;
    end
end
end

totalFastTrackDist = [];
totalSlowTrackDist = [];
for a = 1:length(totalTracks)
if totRotIdx(a) == 1
    if slowTrackIDX(a) > 0
        totalSlowTrackDist = vertcat(totalSlowTrackDist,totalTracks{a}(:,6));
    else
        totalFastTrackDist = vertcat(totalFastTrackDist,totalTracks{a}(:,6));
    end  
end
end

figure;
hold on;
title(['Distance to the cell axis']);
xlabel('Distance (nm)');
ylabel('Probability Density');
edges = 0:10:max(totalDistance);
binSpacing = edges(2)-edges(1);
histogram_counts = histcounts(totalFastTrackDist,(edges+binSpacing/2));
pdf_distribution = histogram_counts;
bar1 = bar(edges(2:end),pdf_distribution,'r');
hold on;
histogram_counts = histcounts(totalSlowTrackDist,(edges+binSpacing/2));
pdf_distribution = histogram_counts;
bar1 = bar(edges(2:end),pdf_distribution,'b');
xlim([0 Lim])
%end

% %plot distance ratio to cell axis for all points
% total_num_dist_ratio = length(totalRatio);
% figure;
% hold on;
% title(['Distance to the cell axis']);
% xlabel('Ratio');
% ylabel('Probability Density');
% edges = 0:0.01:max(totalRatio);
% binSpacing = edges(2)-edges(1);
% histogram_counts = histcounts(totalRatio(clustIDXtot == 0 ).*2,(edges+binSpacing/2));
% pdf_distribution = histogram_counts./(sum(histogram_counts));
% bar1 = bar(edges(2:end),pdf_distribution);
% 
% %plot distance ratio to cell axis for cluster points
% % totalDistanceRatioClust = totalRatio(clustIDXtot == 1);
% totalDistanceRatioClust = totalRatio(clustIDXtot ~= 0);
% total_num_dist_clust = length(totalDistanceRatioClust);
% % figure;
% hold on;
% % title(['Distance to the cell axis']);
% % xlabel('Ratio');
% % ylabel('Probability Density');
% % edges = 0:0.01:max(totalDistanceRatioClust);
% % binSpacing = edges(2)-edges(1);
% histogram_counts = histcounts(totalDistanceRatioClust.*2,(edges+binSpacing/2));
% pdf_distribution = histogram_counts./(sum(histogram_counts));
% bar1 = bar(edges(2:end),pdf_distribution,'r');




%plot cluster center of mass (COM) distance
% total_num_dist = length(totalDistance);
clustDist = zeros(max(clustIDXtot),1);
% for l = 1:max(clustIDX{numCellsTotal})
for l = 1:max(clustIDXtot)
clustDist(l) = mean(totalDistance(clustIDXtot == l));
end
figure;
hold on;
title(['Distance to the cell axis']);
xlabel('Distance (nm)');
ylabel('Probability Density');
edges = 0:10:max(clustDist);
binSpacing = edges(2)-edges(1);
histogram_counts = histcounts(clustDist,(edges+binSpacing/2));
pdf_distribution = histogram_counts./(sum(histogram_counts));
bar1 = bar(edges(2:end),pdf_distribution);

z = 1;









