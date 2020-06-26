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
load v1/rainSTORM_GaussianLS

%% fit checker setup
columns_not_fitted = [1,2,3 ,4,5,6];
%x_column = 2; %what column has x-pos in the return data
%y_column = 3; %what column has y-pos in the return data
n_fits_per_frame = (number_x-size(columns_not_fitted,2))*number_y;

%% fit checker %% clear test
clear res_precision res_accuracy 

% convert all to m

data = SupResParams;
temp = num2cell([SupResParams.x_coord].*mic_pixelsize); % convert to nm
[data.y_coord] = temp{:};
temp = num2cell([SupResParams.y_coord].*mic_pixelsize); % conver to nm
[data.x_coord] = temp{:};

total_fits = 0;

tracked = data;

for i=1:number_x
    for j=1:number_y
        
%     if ismember(i, columns_not_fitted)
%         continue
%     end
    
    temp = data(([data.x_coord] > pos_x(i)-10*mic_pixelsize) & ([data.x_coord] < pos_x(i)+10*mic_pixelsize) & ([data.y_coord] > pos_y(j)-10*mic_pixelsize) & ([data.y_coord] < pos_y(j)+10*mic_pixelsize));
    
    tracked(([tracked.x_coord] > pos_x(i)-10*mic_pixelsize) & ([tracked.x_coord] < pos_x(i)+10*mic_pixelsize) & ([tracked.y_coord] > pos_y(j)-10*mic_pixelsize) & ([tracked.y_coord] < pos_y(j)+10*mic_pixelsize)) =[];
    
    fit_x = [temp.x_coord]';
    fit_y = [temp.y_coord]';
    
    sigma_x = [temp.sig_x]'*mic_pixelsize;
    sigma_y = [temp.sig_y]'*mic_pixelsize;
    
    sigma_x_mean = mean(sigma_x);
    sigma_y_mean = mean(sigma_y);

    sigma_x_std = sum((sigma_x - sigma_x_mean).^2)/(size(sigma_x,1)-1);
    sigma_y_std = sum((sigma_y - sigma_y_mean).^2)/(size(sigma_y,1)-1);

    sigma_mean = mean([sigma_x_mean sigma_y_mean]);
    res_sigma_precision(i,j) = sqrt(sigma_x_std^2 + sigma_y_std^2);
    res_sigma_accuracy(i,j) = sigma_mean - mic_pixelsize;
    
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
    total_fits = total_fits + size(fit_x,1);
    end
end

res_mean_precision = nanmean(res_precision,2);
res_mean_accuracy = nanmean(res_accuracy,2);
res_mean_sigma_precision = nanmean(res_sigma_precision,2);
res_mean_sigma_accuracy = nanmean(res_sigma_accuracy,2);
%% clear test
%clear res_precision res_accuracy res_mean