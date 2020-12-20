clear all
%%
load data_generation_info_v4

SNRs = saved.SNRs;

for i=1:size(SNRs,2)
    
    frame = zeros(9,9);
    gauss_data.I_max=600*SNRs(i);
    gauss_data.sigma=1;
    
    gauss_data.xc = 5;
    gauss_data.yc = 5;
    
    gsize = [9,9];
    [R,C] = ndgrid(1:gsize(1), 1:gsize(2));
    
    u = ((R-gauss_data.xc).^2 + (C-gauss_data.yc).^2)./(2*gauss_data.sigma^2); 
    frame = uint16(gauss_data.I_max*(exp(-u))); 
    
    counts(i) = sum(frame, 'all');

    frame_x = diff(double(frame),1,2);
    frame_x = [frame_x zeros(9,1)]*saved.pixelsize;
    frame_y = diff(double(frame),1,1);
    frame_y = [frame_y; zeros(9,1)'];
    
    frame_double = double(frame);
    
    I = frame_x.*frame_x./frame_double;
    I(isinf(I)|isnan(I)) = 0;
    sigma(i) = 1/sqrt(sum(I,'all'));
end