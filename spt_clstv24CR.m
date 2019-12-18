% spt_clst.m includes several tools to analyze track localization and
% distribution and to present them with types of heat maps. 
%
% Written by Yu-Jen Chen, Department of Anatomy and Structural Biology,
% Albert Einstein College of Medince, Coleman laboratory

% track_summary is at the following format:
% ID; start_Frame; end_Frame; start_X; end_X; start_Y; end_Y; average_X; average_Y 
%
% Updates:
% 10/02/2018 V 2.4 RC added ability to recover centroid of clusters
% 09/21/2015 V 2.3
% Renamed to spt_clst. The output structure is also renamed to 'sptclst'.
% 07/13/2015 V 2.2
% Incorporate 2nd filtering to remove clusters that contain too few tracks
%   (criteria set by mN2)
%
% 07/07/2015 zone_ana 2.1
% Use an octagon instead of a square window for density calculation
% (mdensity4)
%
% binding density heat map with residence time filter
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set parameters for analysis %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = 'Starting row in sptana to start analyzing:   ';
Cell_start = input(prompt); %set the location in sptana structure

prompt = 'Ending row in sptana to complete analysis:   ';
Cell_end = input(prompt); %set the location in sptana structure


mN = 3;             % set the minimum density (contributed by mN number of tracks) for filtering
mN2 = mN;           % sencond filtering with the actual number of tracks withing the high density zone
min_dwt = [0];      % set the minimum dwell time of molecules included in the analysis
enfac = 10;         % set the factor to scale up the map for high resolution (only 10x has been tested)
hdw = 1;            % set hdw for high-resolution density averaging window radius (px)

for cell_id = Cell_start:Cell_end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Import the data for analysis %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tracks = sptana(cell_id).track_summary;
    t_id = sptana(cell_id).nucl_i;
    ndwet = sptana(cell_id).ndwet;
    envlop = sptana(cell_id).envelop;
    acqu = sptana(cell_id).acqu;
    dwet = (tracks(:,3) - tracks(:,2) + 1) * acqu;

    num_inZ = zeros(length(min_dwt), 1);
    trk_inZ = cell(length(min_dwt), 1);
    Z_boundary = cell(length(min_dwt), 1);
    
    for p = 1:length(min_dwt)

         for den_w = hdw
             
            s_i = find(ndwet > min_dwt(p));
            s_tid = t_id(s_i);

            t_X = tracks(s_tid, 8) * enfac;     % scale up!
            t_Y = tracks(s_tid, 9) * enfac;

            %%%%%%%%%%%%%%%%%%%%%%%%
            % Generate density map %
            %%%%%%%%%%%%%%%%%%%%%%%%

            [new_x, new_y] = membound(envlop(:,1), envlop(:,2), enfac);     % scale up!

            mapy = round(max(new_y)+300); mapx = round(max(new_x)+300);
            
            disp(['calculating density map of cell ', num2str(cell_id), '; minimum dwell time: ', num2str(min_dwt(p))]);
            [dmap] = mdensity4(t_X, t_Y, 'nopar', den_w * enfac, new_x, new_y, mapx, mapy);
            devent = dmap * enfac^2 /0.006889/((acqu*2000)/60);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % define and locate clusters %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            m_den = round((mN/(8*(den_w)^2*tan(pi/8)))/0.006889/((acqu*2000)/60), 2); 
                                                        % octagon area is used in normalization

            filweight = devent > m_den;                 % filter the map with minimum density
            hdens = devent .* filweight;                % filtered density map



            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Define clusters through filtered density map %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            dilF = strel('octagon', 3 * floor(den_w * enfac/3));
            dil_filweight = imdilate(filweight, dilF);

            soliZone = bwlabel(dil_filweight);
            periZone = cell(max(max(soliZone)), 1);
            centr = regionprops(soliZone,'Centroid');
            centZone = cat(1, centr.Centroid);

            N_Loc = zeros(1, max(max(soliZone)));
            I_Loc = cell(max(max(soliZone)), 1);

            disp(['assessing clusters of cell ', num2str(cell_id), '; minimum dwell time: ', num2str(min_dwt(p))]);
            
            for i = 1:max(max(soliZone))
                iZone = soliZone == i;

                peri = bwboundaries(iZone, 'noholes');

                periZone{i} = peri{:};
                periZx = periZone{i}(:,1);
                periZy = periZone{i}(:,2);
                [in_id] = inpolygon(t_X, t_Y, periZy, periZx);
                N_Loc(i) = sum(in_id);
                I_Loc{i} = s_tid(in_id);    

            end

            inZ_i = find(N_Loc > mN2);      % remove the high density areas with less than mN2 tracks in them

            num_inZ(p) = length(inZ_i);
            trk_inZ{p} = {I_Loc{inZ_i}};
            Z_boundary{p} = {periZone{inZ_i}};
            Clust_centr{p} = {centr(inZ_i).Centroid};
         end
         
    disp(['calculation completed for cell ', num2str(cell_id), ...
        '; minimum dwell time: ', num2str(min_dwt(p))]);    
    
    end
    
        
    sptclst(cell_id).condition = sptana(cell_id).condition;
    sptclst(cell_id).filename = sptana(cell_id).filename;
    sptclst(cell_id).pathname = sptana(cell_id).pathname;
    
    
    sptclst(cell_id).minimum_dwet = min_dwt;
    sptclst(cell_id).density_window = den_w;
    sptclst(cell_id).num_clusters = num_inZ;
    sptclst(cell_id).boundaries = Z_boundary;
    sptclst(cell_id).clustered_trks = trk_inZ;
    sptclst(cell_id).minimum_den = m_den;
    sptclst(cell_id).minimum_num = mN2;
    sptclst(cell_id).cluster_centroids= Clust_centr;
    
    disp(['data saved for cell ', num2str(cell_id) ]);
end
            
            