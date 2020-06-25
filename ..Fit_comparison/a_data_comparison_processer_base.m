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

pos_x = pos_x - 0.5*mic_pixelsize; % pixel adjust
pos_y = pos_y - 0.5*mic_pixelsize; % pixel adjust

n_frames = 1000;
%% load in data
% load in somehow

columns_not_fitted = [1 2];
x_column = 4; %what column has x-pos in the return data
y_column = 4; %what column has y-pos in the return data
n_fits_per_frame = (number_x-size(columns_not_fitted,2))*number_y;
i_fit = 1; % row index

for i=1:number_x
    for j=1:number_y
    if ismember(i, columns_not_fitted)
        continue
    end
    fit_x = ones(1,n_frames)*59.5+(rand(1,n_frames)-0.5);%data((i-1)*number_y+j:n_fits_per_frame:n_frames*n_fits_per_frame,x_column);
    fit_y = ones(1,n_frames)*19.5+(rand(1,n_frames)-0.5);%data((i-1)*number_y+j:n_fits_per_frame:n_frames*n_fits_per_frame,y_column);
    fit_x_mean = mean(fit_x);
    fit_y_mean = mean(fit_y);
    sigma_x = sum((fit_x - fit_x_mean).^2)/(size(fit_x,2)-1);
    sigma_y = sum((fit_y - fit_y_mean).^2)/(size(fit_y,2)-1);
    res_precision = sqrt(sigma_x^2 + sigma_y^2);
    res_accuracy = sqrt(sum(([pos_x(i) pos_y(j)] - [fit_x_mean fit_y_mean]).^2));
   % dif = 
    end
end