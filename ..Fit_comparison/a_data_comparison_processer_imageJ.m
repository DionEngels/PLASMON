%% setup

load data_generation_info_v4

positions = saved.positions;

mic_pixelsize = saved.pixelsize;

number_x = 20;
number_y = 10;
n_objects = number_x*number_y; %# objects

n_frames = 1000;
%% load in data
% Use "import data"
% Import as Numeric Matrix

%% fit checker setup
x_column = 2; %what column has x-pos in the return data
y_column = 3; %what column has y-pos in the return data

sigma_column = 4;

%% fit checker %% clear test
clear res_precision res_accuracy res_mean res_precision_imageJ
for i=1:size(positions,1)
    pos_x = positions(i,1);
    pos_y = positions(i,2);
    row = floor((i-1)/10)+1;
    y_val = mod(i,number_y);
    if y_val == 0
        y_val = y_val + number_y;
    end
    column = y_val;
    
    temp = data((data(:,x_column) > pos_x-10*mic_pixelsize) & (data(:,x_column) < pos_x+10*mic_pixelsize) & (data(:,y_column) > pos_y-10*mic_pixelsize) & (data(:,y_column) < pos_y+10*mic_pixelsize),:);
    
    fit_x = temp(:,x_column);
    fit_y = temp(:,y_column);
    fit_uncertainty = temp(:,8); % 8 for MLE, 9 for LS

    sigma = temp(:,sigma_column);

    sigma_mean = mean(sigma);

    sigma_std = sum((sigma - sigma_mean).^2)/(size(sigma,1)-1);

    res_sigma_precision(row,column) = sigma_std;
    res_sigma_accuracy(row,column) = sigma_mean - mic_pixelsize;
    
    fit_x_mean = mean(fit_x);
    fit_y_mean = mean(fit_y);
    
    fit_x_std = sum((fit_x - fit_x_mean).^2)/(size(fit_x,1)-1);
    fit_y_std = sum((fit_y - fit_y_mean).^2)/(size(fit_y,1)-1);
    res_precision(row,column) = sqrt(fit_x_std^2 + fit_y_std^2);
    res_accuracy(row,column) = norm([pos_x pos_y] - [fit_x_mean fit_y_mean]);
    res_precision_imageJ(row,column) = mean(fit_uncertainty);
end

res_mean_precision_imageJ = nanmean(res_precision_imageJ,2);
res_mean_precision = nanmean(res_precision,2);
res_mean_accuracy = nanmean(res_accuracy,2);
res_mean_sigma_precision = nanmean(res_sigma_precision,2);
res_mean_sigma_accuracy = nanmean(res_sigma_accuracy,2);
%% clear test
%clear res_precision res_accuracy res_mean