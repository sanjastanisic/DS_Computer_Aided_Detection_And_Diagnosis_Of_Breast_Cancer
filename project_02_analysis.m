clear; clc; close all;
addpath(genpath('./data'));
addpath(genpath('./assets'));
addpath(genpath('./thirdparty-libraries'));

% Create a log file
diary('log.txt');
diary on

% can be '' or '-p01' or '-p02'
sufix = '-p02'; 

dataTable = table;
% add dummy element
dataTable.PatientId = "---";
dataTable.LeftOrRightBrest = "LEFT - RIGHT";
dataTable.ImageView = "CC - MLO";
dataTable.ROIFilePath = "";
dataTable.TestOrTraining = "test - training";
dataTable.BrestDensity = 0;
dataTable.CalcType = "";
dataTable.CalcDistribution = "";
dataTable.Patology = "BENIGN - MALIGNANT";
dataTable.LesionVolume = -1.0;
dataTable.LesionArea = -1.0;
dataTable.SphericalDisproportion = -1.0;
dataTable.Sphericity = -1.0;
dataTable.SurfaceToVolumeRatio = -1.0;

% read test description csv data
filePath = strcat('./data/calc_case_description_test_set', '.csv');
opts = detectImportOptions(filePath);
testSetDesc = readtable(filePath, opts);

% read training description csv data
filePath = strcat('./data/calc_case_description_train_set', '.csv');
opts = detectImportOptions(filePath);
trainSetDesc = readtable(filePath,opts);

%add roi test data
dirPath = "./data/CALC_ROI_Test";
filePath = strcat(dirPath, "/metadata", sufix, ".csv");
opts = detectImportOptions(filePath);
roiTest = readtable(filePath,opts);
for i = 1:size(roiTest,1)
    lf = roiTest{i,16};
    lf = replace( lf, '\', '/');
    lf = replace(lf,'./','/');
    filePath = strcat(dirPath, lf);
    files = dir(strcat(filePath,"/*.dcm"));
    files = struct2table(files);
    if(size(files,1) == 0)
        continue;
    elseif(size(files,1) == 1)
        fileName = files{1, 1};
        filePath = strcat(dirPath, lf, "/", fileName);
    elseif (size(files,1) == 2)
        sby1 = double(files{1, 4});
        sby2 = double(files{2, 4});
        if(sby1>sby2)
            fileName = files{1, 1};
        else
            fileName = files{2, 1};
        end
        filePath = strcat(dirPath, lf, "/", fileName);
    else
        continue;
    end
    try
        info = dicominfo(filePath);
    catch
        continue;
    end
    if isfield(info,'PatientID') == 1
      s = info.PatientID;
      s = replace(s,"Calc-Test_", ""); 
      s = replace(s,".dcm", ""); 
      s = strsplit(s, "_");
      if size(s) < 4
          continue;
      end
      patientId = strcat(s(1), "_", s(2));
      lorb = string(s(3));
      iv = string(s(4));
    else
      continue;
    end  
    %  consult test description csv data
    cond1 = testSetDesc{:,1} == patientId;
    cond2 = testSetDesc{:,3} == lorb;
    cond3 = testSetDesc{:,4} == iv;
    idx = cond1 & cond2 & cond3;
    v = testSetDesc(idx,:);
    if size(v,1) == 0
        continue;
    end
    breastDensity = v{1,2};
    calcType = v{1,7};
    calcDistribution = v{1,8};
    pathology = v{1,10};
    if pathology ~= "BENIGN" &&  pathology ~= "MALIGNANT"
       continue;
    end
    % insert relevant data
    newData = {patientId, lorb, iv, filePath, "test",... 
        breastDensity, calcType, calcDistribution, pathology,...
        -1.0, -1.0, -1.0, -1.0, -1.0};
    dataTable = [dataTable;newData];
end

% add  roi training data
dirPath = "./data/CALC_ROI_Training";
filePath = strcat(dirPath, "/metadata", sufix, ".csv");
opts = detectImportOptions(filePath);
roiTraining = readtable(filePath,opts);
for i = 1:size(roiTraining,1)
    lf = roiTraining{i,16};
    lf = replace( lf, '\', '/');
    lf = replace(lf,'./','/');
    filePath = strcat(dirPath, lf);
    files = dir(strcat(filePath,"/*.dcm"));
    files = struct2table(files);
    if(size(files,1) == 0)
        continue;
    elseif(size(files,1) == 1)
        fileName = files{1, 1};
        filePath = strcat(dirPath, lf, "/", fileName);
    elseif (size(files,1) == 2)
        sby1 = double(files{1, 4});
        sby2 = double(files{2, 4});
        if(sby1>sby2)
            fileName = files{1, 1};
        else
            fileName = files{2, 1};
        end
        filePath = strcat(dirPath, lf, "/", fileName);
    else
        continue;
    end
    try
        info = dicominfo(filePath);
    catch
        continue;
    end
    if isfield(info,'PatientID') == 1
      s = info.PatientID;
      s = replace(s,"Calc-Training_", ""); 
      s = replace(s,".dcm", ""); 
      s = strsplit(s, "_");
      if size(s) < 4
          continue;
      end
      patientId = strcat(s(1), "_", s(2));
      lorb = string(s(3));
      iv = string(s(4));
    else
      continue;
    end
    %  consult training description csv data
    cond1 = trainSetDesc{:,1} == patientId;
    cond2 = trainSetDesc{:,3} == lorb;
    cond3 = trainSetDesc{:,4} == iv;
    idx = cond1 & cond2 & cond3;
    v = trainSetDesc(idx,:);
    if size(v,1) == 0
        continue;
    end
    breastDensity = v{1,2};
    calcType = v{1,7};
    calcDistribution = v{1,8};
    pathology = v{1,10};
    if pathology ~= "BENIGN" &&  pathology ~= "MALIGNANT"
       continue;
    end
    % insert relevant data
    newData = {patientId, lorb, iv, filePath, "training",... 
        breastDensity, calcType, calcDistribution, pathology,...
        -1.0, -1.0, -1.0, -1.0, -1.0};
    dataTable = [dataTable;newData];
end

% delete dummy element from dataTable
dataTable(dataTable.PatientId == "---", :) = [];

% save obtained redults into csv
%writetable(dataTable,strcat('results-f01', sufix, '.csv'));  
% and load it again
%opts = detectImportOptions(strcat('results-f01', sufix, '.csv'));
%dataTable = readtable(strcat('results-f01', sufix, '.csv'), opts);

% process images within dataTable
for i = 1:size(dataTable,1)
    filePath = dataTable{i,4};
    % load image info
    info = dicominfo(filePath);
    if isfield(info,'RescaleSlope') == 1
      slope = info.RescaleSlope;
    else
      slope = 1;
    end
    if isfield(info,'RescaleIntercept') == 1
      intercept = info.RescaleIntercept;
    else
      intercept = 0;
    end
    % calculate the dimension of a single voxel [mm^3]
    ddim = 0.05; % 50 micro meters, from PeerJ paper
    if isfield(info,'PixelSpacing') == 1
       if size(info.PixelSpacing) >= 1
          dimx = info.PixelSpacing(1);
       else
          dimx = ddim; 
       end
       if size(info.PixelSpacing) >= 2
          dimy = info.PixelSpacing(2);
       else
          dimy = ddim;
       end
       if size(info.PixelSpacing) >= 3
          dimz = info.PixelSpacing(3);
       else
          dimz = 1;
       end
    else
       dimx = ddim; 
       dimy = ddim; 
       dimz = 1;
    end
    % load image
    img = dicomread(filePath);
    img = double(img);        
    img = img.*slope + intercept;
    % calculate volume of the leison
    voxel__volume = dimx * dimy * dimz;
    voxel__volume = double(voxel__volume);
    lesion__voxels = sum(img~=0,'all'); % number of voxels with non-zero values
    lesion__volume = lesion__voxels * voxel__volume;
    lesion__volume = double(lesion__volume);
    % calculate area of the lesion
    pixelDimensions = [1 1 1];
    pixelDimensions(1,1) = dimx;
    pixelDimensions(1,2) = dimy;
    pixelDimensions(1,3) = dimz;    
    [lesion__area, surf_mat] = compute__surface(img, double(pixelDimensions));
    lesion__area = double(lesion__area);
    % calculate spherical disproportion
    R_equiv = (lesion__volume*3/(4*pi))^(1/3);
    spherical__disproportion = lesion__area/(4*pi*(R_equiv)^2);
    % calculate sphericity
    sphericity = ((pi^(1/3))*((6*lesion__volume)^(2/3)))/lesion__area;
    % calculate surface-to-volume ratio
    surfacevolume__ratio = lesion__area/lesion__volume;
    dataTable{i,10} = lesion__volume;
    dataTable{i,11} = lesion__area;
    dataTable{i,12} = spherical__disproportion;
    dataTable{i,13} = sphericity;
    dataTable{i,14} = surfacevolume__ratio;
    if mod(i,50) == 0
        writetable(dataTable,strcat('results-f02', sufix, '-', string(i), '.csv'));
    end
end

% save obtained results into csv
writetable(dataTable,strcat('results-f02', sufix, '.csv'));
% and load it again
%opts = detectImportOptions(strcat('results-f02', sufix, '.csv'));
%dataTable = readtable(strcat('results-f02', sufix, '.csv'), opts);

idx = dataTable{:,9} == "BENIGN";
benign__features = dataTable(idx,10:14);
benign__features = table2array(benign__features);
benign__features = double(benign__features);

idx = dataTable{:,9} == "MALIGNANT";
malignant__features = dataTable(idx,10:14);
malignant__features = table2array(malignant__features);
malignant__features = double(malignant__features);

% Unsupervised learning (basic try)
idx = kmeans([benign__features(:,:); malignant__features(:,:)],2,...
    'Replicates',10000);

idx__1_1 = sum(idx(1:size(benign__features,1),:) == 1);
idx__1_2 = sum(idx(1:size(benign__features,1),:) == 2);
idx__2_1 = sum(idx(size(benign__features,1)+1:end,:) == 1);
idx__2_2 = sum(idx(size(benign__features,1)+1:end,:) == 2);

if idx__1_1 < idx__1_2
    temp__a = idx__1_1;
    idx__1_1 = idx__1_2;
    idx__1_2 = temp__a;

    temp__b = idx__2_1;
    idx__2_1 = idx__2_2;
    idx__2_2 = temp__a;
end

unsupervised__accuracy = (idx__1_1 + idx__2_2) / (size(benign__features,1) + size(malignant__features,1));

finito=1;