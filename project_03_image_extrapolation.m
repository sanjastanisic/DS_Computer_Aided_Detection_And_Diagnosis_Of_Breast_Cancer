clear; clc; close all;
addpath(genpath('./data'));
addpath(genpath('./assets'));
addpath(genpath('./thirdparty-libraries'));

% Create a log file
diary('log.txt');
diary on

sufix = '-p03-f01';

opts = detectImportOptions(strcat('results', sufix, '.csv'));
dataTable = readtable(strcat('results', sufix, '.csv'), opts);
 
% process images within dataTable
for i = 1:size(dataTable,1)
    % obtain path of full image
    filePath = dataTable{i,4};
    % load info about full image
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
    % obtain path of mask
    filePathM = dataTable{i,5};
    % load info about mask
    infoM = dicominfo(filePathM);
    if isfield(infoM,'RescaleSlope') == 1
      slopeM = infoM.RescaleSlope;
    else
      slopeM = 1;
    end
    if isfield(infoM,'RescaleIntercept') == 1
      interceptM = infoM.RescaleIntercept;
    else
      interceptM = 0;
    end
    % calculate the dimension of a single voxel [mm^3]
    ddimM = 0.05; % 50 micro meters, from PeerJ paper
    if isfield(infoM,'PixelSpacing') == 1
       if size(infoM.PixelSpacing) >= 1
          dimxM = infoM.PixelSpacing(1);
       else
          dimxM = ddimM; 
       end
       if size(infoM.PixelSpacing) >= 2
          dimyM = infoM.PixelSpacing(2);
       else
          dimyM = ddimM;
       end
       if size(infoM.PixelSpacing) >= 3
          dimzM = infoM.PixelSpacing(3);
       else
          dimzM = 1;
       end
    else
       dimxM = ddimM; 
       dimyM = ddimM; 
       dimzM = 1;
    end
    % load full image
    img = dicomread(filePath);
    img = double(img);        
    img = img.*slope + intercept;
    % load mask
    imgM = dicomread(filePathM);
    imgM = double(imgM);        
    imgM = imgM.*slopeM + interceptM;
    % create extrapolated image
    if size(img,1)== size(imgM,1) && size(img,2)== size(imgM,2) 
      for pi = 1:size(img,1)
        for pj = 1:size(img,2)
            if imgM(pi,pj)==0
                img(pi,pj)=0;
            end
        end
      end
    end
    % save extrapolated image
    newFilePath = replace(filePathM, ".dcm", "-fin.dcm");
    dicomwrite(img,newFilePath); 
end

% save obtained results into csv
writetable(dataTable,strcat('results', sufix, '.csv'));
% and load it again
%opts = detectImportOptions(strcat('results', sufix, '.csv'));
%dataTable = readtable(strcat('results', sufix, '.csv'), opts);

finito=1;