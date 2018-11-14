%%% Charlie Jeynes  %%%%
%%% this file loads in parafin and electronically normalised data, then
%%% groups the data according to whether it is cancerous etc
%%% and then runs a PCA on the average spectra . We were also playing with
%%% the HCA function but this didn't give as good results as the PCA



clc
close
clear all
%%
%filenames = dir('*.mat');  
%% load in data from .mat file
%data = '/Volumes/Minerva_II/Charlie/NM_2CBMA_mucous_all.mat'; 
% data = '/Volumes/Charlie_mac_backup/Charlie Krupaker data oct 2018/NM_2CBMA_mucous_all.mat';
data = '/Volumes/Charlie_mac_backup/Charlie Krupaker data oct 2018/NM_2CBMA_epi.mat';
data = load(data); 

%% %% average all spectra (~20k) to the corresponding ~27 sample names
matrixSampleNum = repmat(data.SampleNum, 1, 801); % 
meanSpectra = []; 
for k = 1:length(data.SAMPLE_LIST)
    averageSpectra = data.DATA(matrixSampleNum == k); 
    lgth = length(matrixSampleNum(matrixSampleNum == k));
    divideby = lgth/801; 
    sampleNumdata = reshape(averageSpectra, divideby, 801); 
    meanSpectra(k, :) = mean(sampleNumdata);
end 
%% PLOT ALL THE SPECTRA JUMBLED COLORS
figure, plot(1:801, meanSpectra);
% legend(data.SAMPLE_LIST{: , 1}) ; 
% figure, gscatter
%% USED THIS BIT TO PLOT THE SPECTRA GROUP BY COLOUR 
% YOU HAVE TO ADJUST THIS BIT AS IS SET TO MUCIN_ALL AT THE MOMENT
T = meanSpectra; 
TMAcancer = T(1:10, :); % HERE SPECIFY THE ROWS FOR EACH GROUP
TMAnormal = T(11:12, :);
normal = T(13:21,:);
tumour = T(22:27,:);
figure, 
hold on
plot(1:801, TMAcancer, 'r', 'DisplayName','TMAcancer')
plot(1:801, TMAnormal, 'g', 'DisplayName','TMAnormal')
plot(1:801, normal, 'b', 'DisplayName','normal' )
plot(1:801, tumour, 'y', 'DisplayName','cancer' )
plot(1:801, T(9, :), 1:801, T(1, :), 1:801, T(7, :), 1:801, T(3, :), 1:801, T(6, :))
legend('show')
hold off

%% PCA on meanSpectra COMBINES MUCIN & AMIDE PEAKS (this can be adjusted)
mucin = meanSpectra(:, 1:300); 
amide = meanSpectra(:, 500:700); 
mucinAmide = [mucin amide]; 
[PCALoadings,score,latent] = pca(mucinAmide);

%% THIS PLOTS PCA 2D FOR PC1 AND PC2 (ALSO KNOWN AS SCORE)
x = score(:, 1); 
y = score(:, 2); 
v = [1 2 3]; % change this to how many groups there are (I looked at data.SAMPLE_LIST)
u = repelem(v,[11 8 6]); % this is how many instances in each group
group = u; 
figure, 
gscatter(x,y,group)
% legend('show')
legend('TM_all','clinNormal', 'clinCancer') %legend('TMcancer', 'TMnormal','clinNormal', 'clinCancer') 



















%% BELOW IS ALL SNIPPETS OF CODE WHICH COULD BE USEFUL
% Y = pdist(meanSpectra); 
% % squareform(Y)
% Z = linkage(Y); 
% figure, dendrogram(Z); 

%% Create a dendrogram 
Y = pdist(meanSpectra,'cityblock');
Z = linkage(Y,'average');
figure
H = dendrogram(Z,'Orientation','left','ColorThreshold',4);
set(H,'LineWidth',2)
c = cophenet(Z,Y)
% yticklabels({'x = 0','x = 5','x = 10'}); 
%xtickslabels(data.SAMPLE_LIST{: , 1}); 
%% Verify 
% c = cophenet(Z,Y)
% %% verify city block method (this is slightly better)
% Y = pdist(meanSpectra,'cityblock');
% Z = linkage(Y,'average');
% c = cophenet(Z,Y)

%% cluster the data
T = cluster(Z,'cutoff',1)
%%
% T = cluster(Z,'maxclust',8)
%length(T)
TMAcancer = T(1:10)
TMAnormal = T(11:12)
normal = T(13:21)
tumour = T(22:27)

% figure, scatter(1:27, T)
%%
figure, 
hold on
scatter(1:27, TMAcancer, 'r', 'DisplayName','TMAcancer')
scatter(1:27, TMAnormal, 'g', 'DisplayName','TMAnormal')
scatter(1:27, normal, 'b', 'DisplayName','normal' )
scatter(1:27, tumour, 'y', 'DisplayName','cancer' )
legend('show')
hold off
%%
% c = clusterdata(meanSpectra,'Linkage','ward','Savememory','on','Maxclust',4);
% grouped = findgroups(alldata, data.SampleNum); 

%% PCA on meanSpectra ALL THE WAVENUMBERS

[coeff,score,latent] = pca(meanSpectra); 

%% PCA on meanSpectra JUST WAVENUMEBR MUCIN
[coeff,score,latent] = pca(meanSpectra(:, 1:200));




%% THIS PLOTS 3 PRINCIPLE COMPOENTS IN 3D

plot(1:401,PCALoadings(:,1:4),'-');
xlabel('Variable');
ylabel('PCA Loading');
legend({'1st Component' '2nd Component' '3rd Component'  ...
	'4th Component'},'location','NW');


figure, 
scatter3(score(:, 1), score(:, 2), score(:, 3)); 
rotate3d on;
%% THIS PLOTS 3 PRINCIPLE COMPOENTS IN 3D WITH GROUPED COLOR

TMAcancer = T(1:10)
TMAnormal = T(11:12)
normal = T(13:21)
tumour = T(22:27)


x = score(:, 1); 
y = score(:, 2); 
z = score(:, 3); 
v = [1 2 3 4];
u = repelem(v,[10 2 9 6])
length(u)
group = u;
colors = brewermap(max(group),'Set1');  %or any other way of creating the colormap
markersize = 20;   %change to suit taste
figure, scatter3(x(:), y(:), z(:), markersize, colors(group,:));
legend('TM', 'tumour', 'tumour'); 
rotate3d on;

%% THIS PLOTS PCA 2D FOR PC1 AND PC2 (ALSO KNOWN AS SCORE)
x = score(:, 1); 
y = score(:, 2); 
v = [1 2 3]; % change this to how many groups there are (I looked at data.SAMPLE_LIST)
u = repelem(v,[11 8 6]); % this is how many instances in each group
group = u; 
figure, 
gscatter(x,y,group)
% legend('show')
legend('TM_all','clinNormal', 'clinCancer') %legend('TMcancer', 'TMnormal','clinNormal', 'clinCancer') 

%% PLOT sample 1 against 16 and tell difference

% PC1_TM_CANCER1 = score(1,:); 
% PC1_clin_norm16 = score(16,:); 

difference = meanSpectra(16,:) - meanSpectra(1,:) ; 

figure, plot(1:801, meanSpectra(1,:), 1:801, meanSpectra(16,:), 1:801, difference); 
legend('show')



%% plot HCA from the 1st 2 components of PCA 

XY = [x y]
Y = pdist(XY,'cityblock');
Z = linkage(Y,'average');
figure
H = dendrogram(Z,'Orientation','left','ColorThreshold',4);
set(H,'LineWidth',2)
c = cophenet(Z,Y)

















%% this is the old krupaker file where H%E sections are lined up by the FTIR image not used here
%% get all the subdirectories
[subfolders] = subdir('/Volumes/Charlie_mac_backup/TAU Project');%subdir('/Users/jcgj201/Documents/MATLAB/FTIR data/Testing_moveResultsTau_data'); % get all the subdirectories
pattern = ["TG" , "WT"]; 
logicalSub = contains(subfolders, pattern); 
subfolders1 = subfolders(logicalSub); 
subfoldersT = subfolders1';

%% load .mat files 
function D = loadData()
D = load('normal_one_EMSC.mat');

end
%%
function [statsTable] = extractData()
 
%% extract data
normData = D.normal_one_EMSC.Data; 
%% squeeze data
normDataSq = squeeze(normData); 
%% plot check
figure, plot(1:801, normDataSq(:, 200)); 
%% transpose for kmeans: rows instances, columns variables
Y = normDataSq'; 
%% Perform k-means
tic
[dataKmean_idx_6,Clustered_spectra] = kmeans(Y,6); 
toc
%% transform kmeans data to image format and show 
x = D.normal_one_EMSC.xy(2:end, 1);
y = D.normal_one_EMSC.xy(2:end, 2);
originalImageSize = zeros(max(y), max(x));
for k = 1 : length(x)
	originalImageSize(y(k), x(k)) = dataKmean_idx_6(k);
end
figure, imagesc(originalImageSize);
cmap = jet(6); 
cmap = flipud(cmap(1:6,:));
% cmap(1,:) = [1,1,1]; set to white, mycolors = [1 0 0; 1 1 0; 0 0 1]; red, green, blue, scale 0,1 
colormap(cmap);
colorbar
axis('on', 'image');
%%
%%%%%% Draw on the kmeans image and get indices%%%%%%%%
%% Put D.data into original image format, so that the ROI getter works
%%%%%% gets inputs for function "ROIcapture" %%%%%%%%
normDataSqT = normDataSq';  %transpose the normalised data so 59000x801 format
x = D.normal_one_EMSC.xy(2:end, 1);
y = D.normal_one_EMSC.xy(2:end, 2);
OriginalData_Image_format = zeros(max(y), max(x), size(normData,1)); % initialise zero array the same size as the original data
for k = 1 : length(x)
	OriginalData_Image_format(y(k), x(k), :) = normDataSqT(k, :);
end
%%
ImgName = D.normal_one_EMSC.FileName; 

%% Call the makeStatsTable on the image
statsTable = makeStatsTable(originalImageSize,OriginalData_Image_format, ImgName); 

end

%% Plot to check that sensible data is coming  - they are all noramlised to 43.9
%%%%%%%%%%%%%this doesn't work for some
%%%%%%%%%%%%%reason%%%%%%%%%%%%%%%%%%%%%%%%%%
% OriginalData_Image_format_Reshape = reshape(OriginalData_Image_format, 65536, 801); 
% figure, plot(1:801, OriginalData_Image_format_Reshape(200, :)); 
% figure, imagesc(OriginalData_Image_format);
% colorbar
%%
function [statsTable] = makeStatsTable(originalImageSize,OriginalData_Image_format, ImgName)
region = []; 
imagename = {}; 
ROIdata = {}; 
statsTable = table(ROIdata, imagename, region); 


% while true
%   drawnow()
%   stop_state = get(YourPushbuttonHandle, 'Value');
%   if stop_state
%     break;
%   end
j = 1;
    while true
    % for k=1:2

        statsTable.ROIdata{j} = ROIcapture(originalImageSize,OriginalData_Image_format); %This functions passes in the image data % name
        statsTable.region(j) = j; 
        statsTable.imagename{j} = ImgName; 
        j = j+1; 
        fig = figure; 
        ax1 = axes('Position',[0 0 1 1],'Visible','off');
        descr = {'click me for another ROI, keyboard press to quit'};
        axes(ax1) % sets ax1 to current axes
        text(.025,0.6,descr)
        w = waitforbuttonpress;
        if w == 0
            disp('Button click')
            disp('continue')
        else
            disp('Key press & break')
            break
        end
    end
end 
%%

function [ROIspectra] = ROIcapture(originalImageSize,OriginalData_Image_format)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is David's impoly bit which allows user to draw a polygon around
% image

%     figure('Name', name); 
%     imagesc(sum(data,3));
    figure, 
    subplot(1,2,2);
    imshow('normal_one_H&E.jpg');
    subplot(1,2,1);
    imagesc(originalImageSize);
    cmap = jet(6); 
    cmap = flipud(cmap(1:6,:));
    % cmap(1,:) = [1,1,1]; set to white, mycolors = [1 0 0; 1 1 0; 0 0 1]; red, green, blue, scale 0,1 
    colormap(cmap);
    colorbar
    axis('on', 'image');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % make it full screen
    
%     ax = axes('Position',[0 0 0.8 1]);
            %"FINISHED-NEXT IMAGE" BUTTON
%         finishButton = uicontrol('String', 'NEXT IMAGE', 'Style', ...
%             'pushbutton', 'Units', 'normalized', 'FontSize', 14);
%         finishButton.Position = [0.8 0.6 0.13 0.05];
%         finishButton.Callback = 'set(gca,''Tag'',''finished'');';
%         finishButton.Enable = 'on';

    ROI  = impoly(gca); %imfreehand imrect(gca);
    wait(ROI); 
    ROIBW = ROI.createMask; 
    [yCoordinates, xCoordinates] = find(ROIBW);

%     saveas(gcf,sprintf('fullImage%d.png',i));
%     close

    x_min = min(xCoordinates); x_max = max(xCoordinates); x_range = x_max-x_min+1;
    y_min = min(yCoordinates); y_max = max(yCoordinates); y_range = y_max-y_min+1;
%     dataROI = -1000*ones(y_range, x_range, size(data,3));
    dataROI = -1000*ones(y_range, x_range, size(OriginalData_Image_format,3)); %intialises a onesROI pixelsxpixelsx801
    for k = 1:numel(yCoordinates)
%         dataROI(yCoordinates(k)-y_min+1, xCoordinates(k)-x_min+1, :) = data(yCoordinates(k), xCoordinates(k), :);
          dataROI(yCoordinates(k)-y_min+1, xCoordinates(k)-x_min+1, :) = OriginalData_Image_format(yCoordinates(k), xCoordinates(k), :);
    end
    dataROI2_im = sum(dataROI,3);
    min_val = min(dataROI2_im(dataROI2_im>=-1000));
    for k1 = 1:y_range
        for k2 = 1:x_range
            if dataROI2_im(k1,k2)<=-1000
                dataROI2_im(k1,k2) = min_val;
            end
        end
    end
    
    ROIspectra = zeros(numel(yCoordinates),size(OriginalData_Image_format,3));
    for k = 1:numel(yCoordinates)
        ROIspectra(k, :) = OriginalData_Image_format(yCoordinates(k), xCoordinates(k), :);
    end

%     ROIspectraMasterNUM1 = size(ROIspectra);
%     ROIspectraMasterNUM2 = ROIspectraMasterNUM1(1); 
%     ROIspectraMasterNUM{i} = ROIspectraMasterNUM2; 
%     ROIspectraMaster{i} = ROIspectra; 
%     ROIspectraMaster1 = cell2mat(ROIspectraMaster'); % create matrix of spectra
%     
%     figure('Name', name);

    figure, imagesc(sum(dataROI2_im,3));
    colorbar
%     saveas(gcf,sprintf('magnified%d.png',i));
end


