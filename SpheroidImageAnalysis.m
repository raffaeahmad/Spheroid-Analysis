% Raffae N. Ahmad
% Virginia Tech
% If using this script please cite:
% "Characterization of Glioma Spheroid Viability and Metastatic Potential 
% Following Monophasic and Biphasic Pulsed Electric Fields,
% Bioelectrochemistry, 2025"

% This code was created to analyze binary masks of spheroids. Two masks should
% be made for the outgrowth region and spheroid body. The file naming
% structure is samplename.tif for original image,
% samplename_mask_outgrowth_dayX.tif for outgrowth mask 
% and samplename_mask_body_dayX.tif for body masks. Output is in pixel
% units and must be converted to um and um^2 if needed. Pixel units can be
% kept if you decide to use relative area between original and
% post-treatment images. 

% Measurements: Number of Protrusions,Area,Centroid, Circularity,
% Eccentricity,MajorAxisLength,MinorAxisLength,MaxFeretProperties,
% MinFeretProperties,Orientation

clc;
clear;
close all;

%% Creates a list of all files that have a tif file extension and are a certain length

%Change to day analyzing ie day0 or day5
day = "day5"; 
experiment = "2023_8_22";

%If set to 1 it will perform max distance math using body mask, set to 0 if
%no body mask exist and only interested in outgrowth metrics
body_calculations = 0;

%If set to 1 will produce a figure for each set of samples overlaying the
%branch points (protrusions) identified from tumor migration boundary and
%save the figure as an image with the same naming convention as original
ovrl_fig = 0;


%Gui select path
folder_path = uigetdir();

% Get a list of all files in the folder
files = dir(fullfile(folder_path, '*.tif'));

% Initialize an empty string array
file_names = string([]);

% Specify the desired length
desired_length = 8; % Change this to your desired length

% Loop through each file and add to the string array if it meets the length and extension criteria
for i = 1:length(files)
    if length(files(i).name) <= desired_length 
        file_names = [file_names; string(files(i).name)];
    end
end

% Display the resulting string array, these are all the files identified
disp('File names with desired length and .tif extension:');
disp(file_names);

%% Using previous list of files opens the BF + masks and finds the centroid

%Image scaling factors from pixels to um, change this based on img prop
scalingfactor = 1.3; %For linear values
sqscalingfactor = scalingfactor^2; %For Area


%All region props of interest
Props = ["Area","Centroid", "Circularity","Eccentricity","MajorAxisLength","MinorAxisLength","MaxFeretProperties","MinFeretProperties","Orientation"];


outgrowthpropstats = table();
bodypropsstats = table();

bcentroid = [];
ocentroid = [];
k = 1;
counter = 1;

for j = 1:length(file_names)

    
    filename = erase(file_names(j),".tif");
    MaxDistances(j,1) = filename;

    %Creates the filenames for the body and outgrowth masks
    bodyfile = filename + "_mask_body_" + day; 
    outgrowthfile = filename + "_mask_outgrowth_" + day;
    
    DeltaAreaMat(j) = 0;

    
    if  ~exist(outgrowthfile + ".tif" , 'file') 
        continue; % Skip to bottom of loop.
    end
   
    Outgrowth = imread(outgrowthfile,"tif");
    outgrowthprops = regionprops(Outgrowth,Props);
    %Finds the xy coordinates of the body centroid 
    
    %%%%%%%%%%%%% Number of protrusions or Branch Points %%%%%%%%%%%%%%%%%%

    centroidsoutgrowth = cat(1,outgrowthprops.Centroid);

    %Idea is calculating the distance between centroid and boundary
    %Coordinates of the outgrowth boundary
    outgrowth_boundary = bwboundaries(Outgrowth);
    
    og_b = cell2mat(outgrowth_boundary(1,1));
            
    %Distance between centroid of outgrowth and boundary coordinates
    Distance = (  ( og_b(:,1) - centroidsoutgrowth(255,2) ).^2 +  ( og_b(:,2)- centroidsoutgrowth(255,1) ).^2 ).^(1/2);
    AverageDist = Distance./mean(Distance);
    
    % This figure creates a line graph of the changing distance and identifies the
    % peaks if you want to visuallize it
    % figure
    % plot(AverageDist) %This figure shows what findPeaks is
    % identifying
    %x = 1:length(AverageDist);
    %findpeaks(AverageDist,x,'MinPeakHeight',0.01,'MinPeakProminence',0.05)

    [pks,loc] = findpeaks(AverageDist,'MinPeakProminence',0.01,'MinPeakHeight',0.05);
    Peaks = og_b(loc,:);
    num_peaks = length(pks);
    
    
    %Displays masks and overlays boundaries, centroid, and peaks
    if ovrl_figure == 1
   
        figure
        imshow(Outgrowth)
        hold on
        %Draws boundary around outgrowth

        for k = 1:length(outgrowth_boundary)
           boundary = outgrowth_boundary{k};
           plot(boundary(:,2), boundary(:,1), 'red', 'LineWidth', 1)
        end

        %Plots centroid and peak location on boundary
        plot(Peaks(:,2),Peaks(:,1),'b*')
        plot(centroidsoutgrowth(255,1),centroidsoutgrowth(255,2),'b*')

        hold off
        A = gcf;
        protrusionfile = filename + "_mask_outgrowth_" + day + "_protrusions.png";
        exportgraphics(A,protrusionfile,"Resolution",150)
        close all
    end


    if body_calculations == 1
        Body = imread(bodyfile,"tif");
        bodyprops = regionprops(Body,Props);
        
        %For some reason this is the quickest way to convert to double
        centroidsbody = cat(1, bodyprops.Centroid);
        centroidsoutgrowth = cat(1,outgrowthprops.Centroid);
        
        %Converts to logicals
        Body_logical = logical(Body); 
        Outgrowth_logical = logical(Outgrowth);
    
        masks(:,:,1) = Body_logical;
        masks(:,:,2) = Outgrowth_logical;
    
        
        %Angle between Centroid of body and centroid of outgrowth
        bcentroid(counter,1) = centroidsbody(255,1);
        bcentroid(counter,2) = centroidsbody(255,2);
        ocentroid(counter,1) = centroidsoutgrowth(255,1);
        ocentroid(counter,2) = centroidsoutgrowth(255,2);
    
        %Change in area
        DeltaArea = outgrowthprops(255,1).Area - bodyprops(255,1).Area;
    
        %Coordinates of the outgrowth boundary
        outgrowth_boundary = bwboundaries(Outgrowth);
        
        og_b = cell2mat(outgrowth_boundary(1,1));
        
        %Distance between centroid of body and boundary coordinates
        Distance = (  ( og_b(:,1) - centroidsbody(255,2) ).^2 +  ( og_b(:,2)- centroidsbody(255,1) ).^2 ).^(1/2);
        
        [Max,Ind] = max(Distance); %Finds max distance and the index of it
        MaxPos = og_b(Ind,:); % Finds the position on the boundary of the max distance


        %Displays masks and overlays boundaries, centroid, and max dist line
        figure
        imshow(Body)
        hold on

        %Draws boundary around outgrowth
        for k = 1:length(outgrowth_boundary)
           boundary = outgrowth_boundary{k};
           plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1)
        end

        %Plots centroid and line from it to the furthest point
        plot(centroidsbody(:,1),centroidsbody(:,2),'b*')
        line([centroidsbody(255,1),MaxPos(2)], [centroidsbody(255,2),MaxPos(1)]);
        hold off
        title(filename)
        B = gcf;
        overlayfile = filename + "_mask_body_" + day +"overlay";
        exportgraphics(B,overlayfile,"Resolution",150)
        close all
    
    
        %Saves the calculated distance in an array, units are in pixels
        %Filename,max dist,body area, outgrowth area, delta area
        MaxDistances(j,4) = string(Max);
        MaxDistances(j,5) = string(bodyprops(255,1).Area);
        MaxDistances(j,6) = string(outgrowthprops(255,1).Area);
        MaxDistances(j,7) = string(DeltaArea);
        DeltaAreaMat(j,1) = DeltaArea;
       
        bodyprops(255).name = filename;
        bodyprops(255).region = 'body';
        bodyprops(255).day = day;
        bodystats = struct2table(bodyprops(255),"AsArray",true);
        bodypropsstats = [bodypropsstats;bodystats];
    end

    outgrowthprops(255).name = filename;
    outgrowthprops(255).region = 'outgrowth';
    outgrowthprops(255).day = day;
    outgrowthprops(255).protrusion = num_peaks;
    outgrowthstats = struct2table(outgrowthprops(255),"AsArray",true);
    
    %creates a table with all stats
    outgrowthpropstats = [outgrowthpropstats;outgrowthstats];
    

    fclose('all');
    counter = counter +1;
end


% Creates excel files
writetable(outgrowthpropstats,experiment + day + '_outgrowth' +'_spheroidstats_.xlsx');
writetable(bodypropsstats,experiment + day +'_body' +'_spheroidstats_.xlsx');    






                                                                                                                                   

