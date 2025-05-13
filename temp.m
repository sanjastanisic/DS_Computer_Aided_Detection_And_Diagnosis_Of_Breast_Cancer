for n = 1:size(files,1)

    load(fullfile(path,files(n).name));

    % Morphological features
        % Volume of the lesion
        voxel__volume = info.PixelDimensions(1,1)/10 * info.PixelDimensions(1,2)/10 * info.PixelDimensions(1,3)/10;
        voxel__volume = double(voxel__volume);
        lesion__voxels = sum(img~=0,'all'); % Returns the number of voxels with non-zero values
        lesion__volume = lesion__voxels * voxel__volume;
        lesion__volume = double(lesion__volume);

        % Surface (area) of the lesion
        [lesion__area, surf_mat] = compute__surface(img, double(info.PixelDimensions));
        lesion__area = double(lesion__area);

        % Spherical disproportion
        R_equiv = (lesion__volume*3/(4*pi))^(1/3);
        spherical__disproportion = lesion__area/(4*pi*(R_equiv)^2);

        % Sphericity
        sphericity = ((pi^(1/3))*((6*lesion__volume)^(2/3)))/lesion__area;

        % Surface-to-volume ratio
        surfacevolume__ratio = lesion__area/lesion__volume;

        morphological__features = [lesion__volume lesion__area spherical__disproportion sphericity surfacevolume__ratio];

        homogeneous__features = [homogeneous__features; morphological__features];
    
    % 

end

% Heterogeneous features
path = '..\lesions\heterogeneous\mat';
files = dir(fullfile(path,'*.mat'));

for n = 1:size(files,1)
   
    load(fullfile(path,files(n).name));
    
    % Volume of the lesion
    lesion__volume = info.Dimensions(1,1)*info.Dimensions(1,1)/10*info.Dimensions(1,2)/10*info.Dimensions(1,3)/10; % num_pixel*vox_dim
    lesion__volume = double(lesion__volume);
    
    % Surface (area) of the lesion
    [lesion__area, surf_mat] = compute__surface(img, info.Dimensions);
    lesion__area = double(lesion__area);
    
    % Spherical disproportion
    R_equiv = (lesion__volume*3/(4*pi))^(1/3);
    spherical__disproportion = lesion__area/(4*pi*(R_equiv)^2);
    
    % Sphericity
    sphericity = ((pi^(1/3))*((6*lesion__volume)^(2/3)))/lesion__area;
    
    % Surface-to-volume ratio
    surfacevolume__ratio = lesion__area/lesion__volume;
    
    morphological__features = [lesion__volume lesion__area spherical__disproportion sphericity surfacevolume__ratio];

    heterogeneous__features = [heterogeneous__features; morphological__features];
    
end

% Unsupervised learning (basic try)
idx = kmeans([homogeneous__features(:,:); heterogeneous__features(:,:)],2,'Replicates',10000);

idx__1_1 = sum(idx(1:size(homogeneous__features,1),:) == 1);
idx__1_2 = sum(idx(1:size(homogeneous__features,1),:) == 2);
idx__2_1 = sum(idx(size(homogeneous__features,1)+1:end,:) == 1);
idx__2_2 = sum(idx(size(homogeneous__features,1)+1:end,:) == 2);

if idx__1_1 < idx__1_2
    temp__a = idx__1_1;
    idx__1_1 = idx__1_2;
    idx__1_2 = temp__a;

    temp__b = idx__2_1;
    idx__2_1 = idx__2_2;
    idx__2_2 = temp__a;
end

unsupervised__accuracy = (idx__1_1 + idx__2_2) / (size(homogeneous__features,1) + size(heterogeneous__features,1));

% Supervised learning (basic try)
labels = [zeros(size(homogeneous__features,1),1); ones(size(heterogeneous__features,1),1)];
SVMModel = fitcsvm([homogeneous__features(:,:); heterogeneous__features(:,:)],labels,'Standardize',true);
CVSVMModel = crossval(SVMModel,'Holdout',0.15);
supervised__accuracy = 1 - kfoldLoss(CVSVMModel);

% Results
disp(['Unsupervised accuracy = ' num2str(unsupervised__accuracy*100) '%']);
disp(['Supervised accuracy = ' num2str(supervised__accuracy*100) '%']);
