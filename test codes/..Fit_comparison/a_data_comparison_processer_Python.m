clear all
%% setup

load data_generation_info_generation_v3_offset_corrected

positions = saved.positions;

mic_pixelsize = saved.pixelsize;
number_y = 10;

n_frames = 1000;
%% load in data
load v9/Gaussian_tol1

%% fit checker setup
x_column = 3; %what column has x-pos in the return data
y_column = 4; %what column has y-pos in the return data

sigma_check = 1;
if sigma_check == 1
    sigma_x_column = 6;
    sigma_y_column = 7;
end

%% fit checker %% clear test
clear res_precision res_accuracy

data = Localizations;
data(:,x_column) = (data(:, x_column))*mic_pixelsize; % convert to nm, compensate for pixel offset MATLAB
data(:,y_column) = (data(:, y_column))*mic_pixelsize; % convert to nm, compensate for pixel offset MATLAB

total_fits = 0;

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
    
    if sigma_check == 1
        sigma_x = temp(:,sigma_x_column)*mic_pixelsize; % convert to nm
        sigma_y = temp(:,sigma_y_column)*mic_pixelsize; % convert to nm
        
        sigma_x_mean = mean(sigma_x);
        sigma_y_mean = mean(sigma_y);
        
        sigma_x_std = std(sigma_x);
        sigma_y_std = std(sigma_y);
        
        sigma_mean = mean([sigma_x_mean sigma_y_mean]);
        res_sigma_precision(row,column) = sqrt(sigma_x_std^2 + sigma_y_std^2);
        res_sigma_accuracy(row,column) = sigma_mean - mic_pixelsize;
        
        intensity = temp(:, 5);
        mean_intensity = mean(intensity);
        mean_background = mean(temp(:,8));
        counts = mean_intensity*2*pi;
        res_crlb(row, column) = sqrt(sigma_mean^2/counts+mic_pixelsize^2/12/counts + 4*sqrt(pi)*mic_pixelsize^3*mean_background/mic_pixelsize/counts^2);
    end
        
    
    
    
    fit_x_mean = mean(fit_x);
    fit_y_mean = mean(fit_y);
    
    fit_x_std = std(fit_x);
    fit_y_std = std(fit_y);
    res_precision(row,column) = sqrt(fit_x_std^2 + fit_y_std^2);
    res_accuracy(row,column) = norm([pos_x pos_y] - [fit_x_mean fit_y_mean]);
    
        %figure
%         hold on
%         scatter(fit_x, fit_y)
%         scatter(pos_x,  pos_y, 'x','g', 'LineWidth',5)
%         scatter(fit_x_mean,  fit_y_mean, 'x','r', 'LineWidth',5)
%         total_fits = total_fits + size(fit_x,1);
end

res_mean_precision = nanmean(res_precision,2);
res_mean_accuracy = nanmean(res_accuracy,2);
if sigma_check == 1
    res_mean_crlb = nanmean(res_crlb,2);
    res_mean_sigma_precision = nanmean(res_sigma_precision,2);
    res_mean_sigma_accuracy = nanmean(res_sigma_accuracy,2);
end
%% clear test
%clear res_precision res_accuracy res_mean