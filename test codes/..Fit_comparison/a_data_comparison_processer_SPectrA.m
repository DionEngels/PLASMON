clear all
%% setup
load data_generation_info_generation_v3_offset_corrected

positions = saved.positions;

pos_x = positions(:,1);
pos_y = positions(:,2);

mic_pixelsize = saved.pixelsize;
number_x = 20;
number_y = 10;

n_frames = 1000;
%% load in data
load v7_z_dataset_v3/SPectrA

%% fit checker setup
columns_not_fitted = [1];
n_fits_per_frame = (number_x-size(columns_not_fitted,2))*number_y;
i_fit = 1; % row index
%% normal
x_column = 3; %what column has x-pos in the return data
y_column = 5; %what column has y-pos in the return data
sigma_column = 4;
%% log
% x_column = 2; %what column has x-pos in the return data
% y_column = 4; %what column has y-pos in the return data
% sigma_column = 3;
%% fit checker %% clear test
clear res_precision res_accuracy 

row = 0;
column = 0;

data = [];

for i=1:size(ParticleFile,2)
    fit_x = (ParticleFile(i).newField.twoDGauss(:,x_column)+ParticleFile(i).newField.Location(1)-0.5)*mic_pixelsize; % convert to nm
    fit_y = (ParticleFile(i).newField.twoDGauss(:,y_column)+ParticleFile(i).newField.Location(2)-0.5)*mic_pixelsize; % convert to nm
    
    sigma = (ParticleFile(i).newField.twoDGauss(:,sigma_column))*mic_pixelsize; % convert to nm
    
    add = [fit_x fit_y sigma];
    
    data = [data ; add];
    
end

for i=1:size(positions,1)
    pos_x = positions(i,1);
    pos_y = positions(i,2);
    row = floor((i-1)/10)+1;
    y_val = mod(i,number_y);
    if y_val == 0
        y_val = y_val + number_y;
    end
    column = y_val;
    
    temp = data((data(:,1) > pos_x-10*mic_pixelsize) & (data(:,1) < pos_x+10*mic_pixelsize) & (data(:,2) > pos_y-10*mic_pixelsize) & (data(:,2) < pos_y+10*mic_pixelsize),:);
    
    fit_x = temp(:,1);
    fit_y = temp(:,2);

    sigma = temp(:,3);

    sigma_mean = mean(sigma);

    sigma_std = std(sigma);

    res_sigma_precision(row,column) = sigma_std;
    if res_sigma_precision(row,column) > 200
        res_sigma_precision(row,column) = NaN;
    end
    res_sigma_accuracy(row,column) = sigma_mean - mic_pixelsize;
    
    fit_x_mean = mean(fit_x);
    fit_y_mean = mean(fit_y);
    
    fit_x_std = std(fit_x);
    fit_y_std = std(fit_y);
    res_precision(row,column) = sqrt(fit_x_std^2 + fit_y_std^2);
    if res_precision(row,column) > 200
        res_precision(row,column) = NaN;
    end
    res_accuracy(row,column) = norm([pos_x pos_y] - [fit_x_mean fit_y_mean]);
end
%%
res_mean_precision = nanmean(res_precision,2);
res_mean_accuracy = nanmean(res_accuracy,2);
res_mean_sigma_precision = nanmean(res_sigma_precision,2);
res_mean_sigma_accuracy = nanmean(res_sigma_accuracy,2);
%% clear test
%clear res_precision res_accuracy res_mean