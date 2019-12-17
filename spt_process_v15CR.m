% This program is designed reprocess the spt data for further analysis 
% with analytical tools in the section, "analysis". 
% 
% The programe takes the output from evalSPT and carry out the following
% works:
% a. summarize track information (starting and ending coordinates and etc.)
% b. define residence time for each track
% c. manually define nuclear envelop and filter out tracks outside the NE
% d. save all processed information into structure array for further use
%
%User must enter in parameters and preset parameters
%User then loads file from evalSPT
%After Figure 102 is displayed user must define the nucleus using 25 mouse
%clicks approximating the nuclear envelop.
%Information is then saved in sptana structure at position defined by
%cell_id set in parameters.

%%% Updates:
% 7/7/2015 15b
% Renamed to spt_process v.15
% Debug after updates in functions, membound and inbound.
% Set the magnifying factor for both membound and inbound to get high
% resolution map in following analysis.
%
% 5/29/2015 14b
% Clened up and separated some tools, such as stage shift evaluation. 
%
% 5/14/2015 14b  
% Updated the section that plots all molecules with indicated dwell time.
% The dwell time is now color coded. 
% Density for heat map is now calculated in a custom function 'mdensity'. 
% 2D molecule displacement is calculated by the function 'dmove2d' for each
% track.
%
%clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%
% set the parameters %
%%%%%%%%%%%%%%%%%%%%%%
prompt = 'Please Define the condition and Cell Identifier:   ';
Cell_Identifier = input(prompt, 's'); %set the location in sptana structure

prompt = 'Please Define a number (row) where to put data in sptana structure:   ';
cell_id = input(prompt); %set the location in sptana structure

prompt = 'Please Define the total number of frames in movie:   ';
total_frame = input(prompt); % set the total number of frames in the movie

prompt = 'Please Define the acquisition time in seconds:   ';
acqu = input(prompt); % set the acquistion time + mechanical delay in seconds

prompt = 'Please Define the rate of photobleaching:   ';
photobleach = input(prompt); % set the acquistion time + mechanical delay in seconds
plotfig = 1;

%%%%%%%%%%%%%%%%%%%%%
% preset parameters %
%%%%%%%%%%%%%%%%%%%%%

p = cell_id;      % this is to set the position in the structure array for saving the current data

sptana(p).condition = Cell_Identifier;              % set manually
sptana(p).condition = sptana(p).condition;
sptana(p).diffusion = 0.05;                             % set manually from slimFAST parameters used
sptana(p).off_frame = 3;                                %  set manually from slimFAST parameters used
sptana(p).slimfast_window = 11;                          %  set manually from slimFAST parameters used
sptana(p).bleachrate = photobleach;                              % set manually from photobleachV1 program

sptana(p).total_frame = total_frame;                    % set manually
sptana(p).acqu = acqu;                                  % set manually

%%%%%%%%%%%%%%%%%%%
% import spt data %
%%%%%%%%%%%%%%%%%%%

[FileName, Pathname] = uigetfile('.txt');
file=[Pathname, FileName];
data = importdata(file);
ind_hits = unique(data(:,4));   
N_hits = length(ind_hits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build the spt summary & define residence time %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tracklength = zeros(N_hits,1); 
track_start_f = zeros(N_hits,1); track_end_f = zeros(N_hits,1);
track_start_x = zeros(N_hits,1); track_start_y = zeros(N_hits,1);
track_end_x = zeros(N_hits,1); track_end_y = zeros(N_hits,1);
avgd = zeros(N_hits,2);
tracks = zeros(N_hits,9);

for i=1:N_hits
    logi=find(data(:,4)==ind_hits(i));
    tracklength(i) = data(logi(end),3)-data(logi(1),3)+1;
    track_start_f(i) = data(logi(1),3); track_end_f(i) = data(logi(end),3);
    track_start_x(i) = data(logi(1),1); track_start_y(i) = data(logi(1),2);
    track_end_x(i) = data(logi(end),1); track_end_y(i) = data(logi(end),2);
    avgd(i, 1) = mean(data(logi, 1)); avgd(i, 2) = mean(data(logi, 2)); 
    tracks(i,:) = [i,track_start_f(i), track_end_f(i), track_start_x(i), track_end_x(i),...
        track_start_y(i), track_end_y(i), avgd(i, 1), avgd(i, 2)];    
end

dwet = acqu .* tracklength;         % residence time defined

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quick look at the spt distribution %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotfig == 1
    
    dt_range = [10 50];

    figure(101); clf;
    plot(avgd(:,1), avgd(:,2), 'rx', 'LineWidth', 2.5, 'MarkerSize', 8);
    title(['N spot = ', num2str(N_hits)]);
    axis equal; axis tight;
    set(gca,'YDir','reverse');
    setxlim = get(gca, 'xlim');
    setylim = get(gca, 'ylim');

    dt_get = (dwet>dt_range(1)) & (dwet<dt_range(2));
    figure(102); clf;
    plot(avgd(:,1), avgd(:,2), 'b.', 'MarkerSize', 12); hold on;
    plot(avgd(dt_get, 1), avgd(dt_get, 2), 'rx', 'LineWidth', 2.5, 'MarkerSize', 8); hold off;
    axis equal; axis tight; set(gca,'YDir','reverse');
    set(gca, 'xlim', setxlim, 'ylim', setylim);
    
end

%%  Tracing out nuclear membrane and determine spatial distribution of dwell times

figure(102)
[mem_x, mem_y] = ginput(25);        % repeat this until the NE is satisfactory

% mem_x = sptana1(cell_id).envelop(:,1);
% mem_y = sptana1(cell_id).envelop(:,2);
setxlim = [min(mem_x)-5, max(mem_x)+5];
setylim = [min(mem_y)-5, max(mem_y)+5];


%% create the continuous envelop
% Run this section first after the defining points have been set.
% This envelope is required for all the following analysis.
enfac = 1;      % magnify the image size for higher resolution

[new_x, new_y] = membound(mem_x, mem_y, enfac); 

figure(21);clf;

scatter(mem_x * enfac, mem_y * enfac, 24 ,'bx'); hold on;
scatter(new_x, new_y, 4, 'b', 'fill'); hold off;

axis equal; axis tight;
set(gca, 'Ydir', 'reverse');
set(gca, 'xlim', setxlim * enfac, 'ylim', setylim * enfac);


%% verify molecules locating inside the nuclear envelope
%
enfac = 10;

xc = avgd(:,1);     % Set array for the original identified molecules.
yc = avgd(:,2);
be_id = ind_hits;
% nMol_x = [];        % X (or column) coordinates for molecules verified being in the nucleus.
% nMol_y = [];        % Y (or row) coordinates ......
% ndwet = [];         % molecules within the defined membrane 
% nucl_i = [];        % index info for molecules inside the nucleus (index to dwet matrix or avgd)

[nMol_x, nMol_y ,nucl_i] = inbound4(xc, yc, be_id, new_x, new_y, enfac);
ndwet = dwet(nucl_i);

[Env_i, nucA] = boundary(nMol_x, nMol_y);
[Env_x, Env_y] = membound(nMol_x(Env_i), nMol_y(Env_i), 1);

figure(201); clf;
plot(nMol_x, nMol_y, 'r.', 'LineWidth', 2, 'MarkerSize',10); hold on;
scatter(Env_x, Env_y, 8 ,'b','fill'); hold off;
title(['Number of tracks: ', num2str(length(nucl_i))]);
axis equal; axis tight;
set(gca,'YDir','reverse');
set(gca, 'xlim', setxlim, 'ylim', setylim);

% rebuild track data with nuclear tracks only

nucl_data = cell(length(nucl_i), 1);

for i = 1:length(nucl_i)
    n_i = find(data(:,4) == nucl_i(i));
    nucl_data{i} = data(n_i, :);
end

disp('nuclear data filtered');


%% save the processed track data in the structure array
%
% Update sptana from the old version
p = cell_id;


nucl_data = cell(length(nucl_i), 1);

for i = 1:length(nucl_i)
    n_i = find(data(:,4) == nucl_i(i));
    nucl_data{i} = data(n_i, :);
end

sptana(p).nucl_data = nucl_data;
sptana(p).nucl_summary = tracks(nucl_i, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sptana(p).filename = FileName;
sptana(p).pathname = Pathname;
sptana(p).envelop = [mem_x, mem_y];

sptana(p).data = data;
sptana(p).num_track = N_hits;
sptana(p).track_summary = tracks;

sptana(p).nucl_data = nucl_data;
sptana(p).nucl_i = nucl_i;
sptana(p).ndwet = ndwet;
sptana(p).nucl_summary = tracks(nucl_i, :);

disp('data saved');
%% save the processed data in file

uisave({'sptana'});
clearvars -except sptana
