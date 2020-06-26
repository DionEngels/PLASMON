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

pos_x = (pos_x - 0.5*mic_pixelsize); % pixel adjust, in nm
pos_y = (pos_y - 0.5*mic_pixelsize); % pixel adjust, in nm

n_frames = 1000;
%% load in data
load v2/Python_ScipyPhasor_Localizations

%% fit checker setup
x_column = 4; %what column has x-pos in the return data
y_column = 3; %what column has y-pos in the return data

sigma_check = 1;
if sigma_check == 1
    sigma_x_column = 6;
    sigma_y_column = 7;
end

%% fit checker %% clear test
clear res_precision res_accuracy

data = Localizations;
data(:,x_column) = data(:, x_column)*mic_pixelsize; % convert to nm
data(:,y_column) = data(:, y_column)*mic_pixelsize; % convert to nm
data = data(data(:,x_column)>0,:);

total_fits = 0;

for i=1:number_x
    for j=1:number_y
        
        temp = data((data(:,x_column) > pos_x(i)-10*mic_pixelsize) & (data(:,x_column) < pos_x(i)+10*mic_pixelsize) & (data(:,y_column) > pos_y(j)-10*mic_pixelsize) & (data(:,y_column) < pos_y(j)+10*mic_pixelsize),:);
        
        
        fit_x = temp(:,x_column);
        fit_y = temp(:,y_column);
        
        if sigma_check == 1
            sigma_x = temp(:,sigma_x_column)*mic_pixelsize; % convert to nm
            sigma_y = temp(:,sigma_y_column)*mic_pixelsize; % convert to nm

            sigma_x_mean = mean(sigma_x);
            sigma_y_mean = mean(sigma_y);
            
            sigma_x_std = sum((sigma_x - sigma_x_mean).^2)/(size(sigma_x,1)-1);
            sigma_y_std = sum((sigma_y - sigma_y_mean).^2)/(size(sigma_y,1)-1);
            
            sigma_mean = mean([sigma_x_mean sigma_y_mean]);
            res_sigma_precision(i,j) = sqrt(sigma_x_std^2 + sigma_y_std^2);
            res_sigma_accuracy(i,j) = sigma_mean - mic_pixelsize;
        end
        
        
        
        
        fit_x_mean = mean(fit_x);
        fit_y_mean = mean(fit_y);
        
        fit_x_std = sum((fit_x - fit_x_mean).^2)/(size(fit_x,1)-1);
        fit_y_std = sum((fit_y - fit_y_mean).^2)/(size(fit_y,1)-1);
        res_precision(i,j) = sqrt(fit_x_std^2 + fit_y_std^2);
        res_accuracy(i,j) = norm([pos_x(i) pos_y(j)] - [fit_x_mean fit_y_mean]);
        
        %     figure
        %     scatter(fit_x, fit_y)
        %     hold on
        %     scatter(pos_x(i),  pos_y(j), 'x','r', 'LineWidth',5)
        %     total_fits = total_fits + size(fit_x,1);
    end
end

res_mean_precision = nanmean(res_precision,2);
res_mean_accuracy = nanmean(res_accuracy,2);
res_mean_sigma_precision = nanmean(res_sigma_precision,2);
res_mean_sigma_accuracy = nanmean(res_sigma_accuracy,2);
%% clear test
%clear res_precision res_accuracy res_mean