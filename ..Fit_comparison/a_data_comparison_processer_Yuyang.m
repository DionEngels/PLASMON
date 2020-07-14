%clear all
close all
%% load in data
load v7_Yuyang/Last_Fit_7x7
load v7_Yuyang/metadata
load v7_Yuyang/ROI_locations


%% fit checker setup
x_column = 4; %what column has x-pos in the return data
y_column = 3; %what column has y-pos in the return data
x_column_Thunder = 2;
y_column_Thunder = 3;

%% fit checker %% clear test

data = Localizations;
% data(:,x_column) = (data(:, x_column))*mic_pixelsize; % convert to nm, compensate for pixel offset MATLAB
% data(:,y_column) = (data(:, y_column))*mic_pixelsize; % convert to nm, compensate for pixel offset MATLAB
%data = data(data(:,x_column)<0 | data(:,x_column)> height ,:);

data_ThunderSTORM = ThunderSTORMresults;
% data_ThunderSTORM(:,x_column_Thunder) = (data_ThunderSTORM(:, x_column_Thunder))/2; % convert to nm, compensate for pixel offset MATLAB
% data_ThunderSTORM(:,y_column_Thunder) = (data_ThunderSTORM(:, y_column_Thunder))/2; % convert to nm, compensate for pixel offset MATLAB


x = 374 - x;
% figure
% scatter(y,x)
%%
figure
hold on
scatter(data(:,x_column), data(:,y_column))
%scatter(data_ThunderSTORM(:,x_column_Thunder), data_ThunderSTORM(:,y_column_Thunder))
%%
hold off
figure
frames = 1:sequence_count;
x_pos = data(data(:,2)==123,x_column);
y_pos = data(data(:,2)==123,y_column);


intensity = data(data(:,2)==124,5);
plot(frames,intensity)


% 
% total_fits = 0;
% 
% for i=1:size(positions,1)
%     
%     pos_x = positions(i,1);
%     pos_y = positions(i,2);
%     row = floor((i-1)/10)+1;
%     y_val = mod(i,number_y);
%     if y_val == 0
%         y_val = y_val + number_y;
%     end
%     column = y_val;
%     
%     temp = data((data(:,x_column) > pos_x-10*mic_pixelsize) & (data(:,x_column) < pos_x+10*mic_pixelsize) & (data(:,y_column) > pos_y-10*mic_pixelsize) & (data(:,y_column) < pos_y+10*mic_pixelsize),:);
%     
%     fit_x = temp(:,x_column);
%     fit_y = temp(:,y_column);
%     
%     if sigma_check == 1
%         sigma_x = temp(:,sigma_x_column)*mic_pixelsize; % convert to nm
%         sigma_y = temp(:,sigma_y_column)*mic_pixelsize; % convert to nm
%         
%         sigma_x_mean = mean(sigma_x);
%         sigma_y_mean = mean(sigma_y);
%         
%         sigma_x_std = std(sigma_x);
%         sigma_y_std = std(sigma_y);
%         
%         sigma_mean = mean([sigma_x_mean sigma_y_mean]);
%         res_sigma_precision(row,column) = sqrt(sigma_x_std^2 + sigma_y_std^2);
%         res_sigma_accuracy(row,column) = sigma_mean - mic_pixelsize;
%     end
%     
%     
%     
%     
%     fit_x_mean = mean(fit_x);
%     fit_y_mean = mean(fit_y);
%     
%     fit_x_std = std(fit_x);
%     fit_y_std = std(fit_y);
%     res_precision(row,column) = sqrt(fit_x_std^2 + fit_y_std^2);
%     res_accuracy(row,column) = norm([pos_x pos_y] - [fit_x_mean fit_y_mean]);
%     
%         %figure
%         hold on
%         scatter(fit_x, fit_y)
%         scatter(pos_x,  pos_y, 'x','g', 'LineWidth',5)
%         scatter(fit_x_mean,  fit_y_mean, 'x','r', 'LineWidth',5)
%         total_fits = total_fits + size(fit_x,1);
% end
% 
% res_mean_precision = nanmean(res_precision,2);
% res_mean_accuracy = nanmean(res_accuracy,2);
% if sigma_check == 1
%     res_mean_sigma_precision = nanmean(res_sigma_precision,2);
%     res_mean_sigma_accuracy = nanmean(res_sigma_accuracy,2);
% end
% %% clear test
% %clear res_precision res_accuracy res_mean