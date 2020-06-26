clear all
%% setup
number_x = 20;
number_y = 10;
n_objects = number_x*number_y; %# objects

SNRs = [0.25, 0.5,  0.75, 1:1:10, 15, 20 30 40 60 80 100];

pixel_spacing_x = 20;
pixel_spacing_y = pixel_spacing_x;

mic_pixels_x = (number_x+1)*pixel_spacing_x; %# pixels
mic_pixels_y = (number_y+1)*pixel_spacing_y; %# pixels

mic_pixelsize = 200;

pos_x = [pixel_spacing_x:pixel_spacing_x:number_x*pixel_spacing_x]*mic_pixelsize;
pos_y = [pixel_spacing_y:pixel_spacing_y:number_y*pixel_spacing_y]*mic_pixelsize;

pos_x = (pos_x - 0.5*mic_pixelsize); % pixel adjust
pos_y = (pos_y - 0.5*mic_pixelsize); % pixel adjust

n_frames = 1000;
%% load in data
load v1/SPectrA_GaussianLS

%% fit checker setup
columns_not_fitted = [1];
x_column = 3; %what column has x-pos in the return data
y_column = 5; %what column has y-pos in the return data
n_fits_per_frame = (number_x-size(columns_not_fitted,2))*number_y;
i_fit = 1; % row index

sigma_column = 4;

%% fit checker %% clear test
clear res_precision res_accuracy 

row = 0;
column = 0;

for i=1:(number_x-size(columns_not_fitted,1))*number_y
    
    if mod(i-1,number_x-size(columns_not_fitted,1)) == 0
        row = row + 1;
        column = size(columns_not_fitted,1);
    end
    column =  column + 1;
    
    %close all
    fit_x = (ParticleFile(i).newField.twoDGauss(:,x_column)+ParticleFile(i).newField.Location(1)-0.5)*mic_pixelsize; % convert to nm
    fit_y = (ParticleFile(i).newField.twoDGauss(:,y_column)+ParticleFile(i).newField.Location(2)-0.5)*mic_pixelsize; % convert to nm
    
    sigma = (ParticleFile(i).newField.twoDGauss(:,sigma_column))*mic_pixelsize; % convert to nm
    
    sigma_mean = mean(sigma);

    res_sigma_precision(column,row) = sum((sigma - sigma_mean).^2)/(size(sigma,1)-1);
    res_sigma_accuracy(column,row) = sigma_mean - mic_pixelsize;

    fit_x_mean = mean(fit_x);
    fit_y_mean = mean(fit_y);
    fit_x_std = sum((fit_x - fit_x_mean).^2)/(size(fit_x,1)-1);
    fit_y_std = sum((fit_y - fit_y_mean).^2)/(size(fit_y,1)-1);
    
    res_precision(column,row) = sqrt(fit_x_std^2 + fit_y_std^2);
    res_accuracy(column,row) = sqrt(sum(([pos_x(column) pos_y(row)] - [fit_x_mean fit_y_mean]).^2));
    if column == 19
        continue
    end
    %figure
    %scatter(fit_x, fit_y)
    %hold on
    %scatter(pos_x(column),  pos_y(row), 'x','r', 'LineWidth',5)
end

res_mean_precision = nanmean(res_precision,2);
res_mean_accuracy = nanmean(res_accuracy,2);
res_mean_sigma_precision = nanmean(res_sigma_precision,2);
res_mean_sigma_accuracy = nanmean(res_sigma_accuracy,2);
%% clear test
%clear res_precision res_accuracy res_mean