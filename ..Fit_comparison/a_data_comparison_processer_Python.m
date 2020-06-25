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

pos_x = (pos_x - 0.5*mic_pixelsize)/1e9; % pixel adjust
pos_y = (pos_y - 0.5*mic_pixelsize)/1e9; % pixel adjust

n_frames = 1000;
%% load in data
load v2/Python_ScipyPhasor_highest_tol_Localizations

%% fit checker setup
x_column = 4; %what column has x-pos in the return data
y_column = 3; %what column has y-pos in the return data

%% fit checker %% clear test
clear res_precision res_accuracy 

data = Localizations;
data(:,x_column) = data(:, x_column)*mic_pixelsize/1e9;
data(:,y_column) = data(:, y_column)*mic_pixelsize/1e9;
data = data(data(:,x_column)>0,:);

total_fits = 0;

for i=1:number_x
    for j=1:number_y
        
%     if ismember(i, columns_not_fitted)
%         continue
%     end
    
    temp = data((data(:,x_column) > pos_x(i)-10*mic_pixelsize/1e9) & (data(:,x_column) < pos_x(i)+10*mic_pixelsize/1e9) & (data(:,y_column) > pos_y(j)-10*mic_pixelsize/1e9) & (data(:,y_column) < pos_y(j)+10*mic_pixelsize/1e9),:);
    
    %tracked(([tracked.x_coord] > pos_x(i)-10*mic_pixelsize/1e9) & ([tracked.x_coord] < pos_x(i)+10*mic_pixelsize/1e9) & ([tracked.y_coord] > pos_y(j)-10*mic_pixelsize/1e9) & ([tracked.y_coord] < pos_y(j)+10*mic_pixelsize/1e9)) =[];
    
    fit_x = temp(:,x_column);
    fit_y = temp(:,y_column);
    
    fit_x_mean = mean(fit_x);
    fit_y_mean = mean(fit_y);
    sigma_x = sum((fit_x - fit_x_mean).^2)/(size(fit_x,1)-1);
    sigma_y = sum((fit_y - fit_y_mean).^2)/(size(fit_y,1)-1);
    res_precision(i,j) = sqrt(sigma_x^2 + sigma_y^2);
    res_accuracy(i,j) = sqrt(sum(([pos_x(i) pos_y(j)] - [fit_x_mean fit_y_mean]).^2));
    
%     figure
%     scatter(fit_x, fit_y)
%     hold on
%     scatter(pos_x(i),  pos_y(j), 'x','r', 'LineWidth',5)
%     total_fits = total_fits + size(fit_x,1);
    end
end

% convert all back to nm

res_precision = res_precision*1e9;
res_accuracy = res_accuracy*1e9;

res_precision_mean = nanmean(res_precision,2);
res_accuracy_mean = nanmean(res_accuracy,2);
%% clear test
%clear res_precision res_accuracy res_mean