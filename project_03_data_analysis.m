clear; clc; close all;
addpath(genpath('./data'));
addpath(genpath('./assets'));
addpath(genpath('./thirdparty-libraries'));

% Create a log file
diary('log.txt');
diary on

sufix = '-p03-part'; % '-p01-f01'

opts = detectImportOptions(strcat('results', sufix, '.csv'));
dataTable = readtable(strcat('results', sufix, '.csv'), opts);

idx = dataTable{:,10} == "BENIGN";
benign__features = dataTable(idx,11:15);

idx = dataTable{:,10} == "MALIGNANT";
malignant__features = dataTable(idx,11:15);

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