function [Population] = LocateBeadsProfile(ImageFrame, BeadImage, Diagnostic, MinimumCorrelationFraction, MaximumLocalCorrelationFraction, EdgeWidth)
% Function LocateBeadsProfile
%   Image Analysis Framework 
%   Version 1.6 (Emiel Visser)
% This function will find the particles (beads) in an image ImageFrame
% using an example particle BeadImage.
%
% Typical usage: 
%         LocateBeadsProfile(FrameData, BeadFilter, false, 0.005, 0.75);
%
% * Input parameters:
%         ImageFrame ::    The image to be filtered (MxN matrix)
%         BeadImage ::     The particle image used to find the particles
%         Diagnostic ::    Default: false. If set to true you'll get a UI to
%                            play with the thresholds
%         MinimumCorrelationFraction ::     
%                          Value between 0 and 1. This threshold sets the
%                            minimum intensity a particle is required to 
%                            have to be detected. This is relative to the
%                            maximum observed correlation between the
%                            imageframe and the beadimage.
%         MaximumLocalCorrelationFraction ::
%                          Value between 0 and 1. Rejects 'weirdly' shaped
%                            particle. A check is done on the
%                            local environment of the found particle
%                            locations. If a correlation is found above 
%                            MaximumLocalCorrelationFraction * the
%                            correlation value of the found particle, it is
%                            discarded. See 'LocalMaximumMask_Ring' in code
% * Output parameters: 
%         Population ::     Structure containing the number of found
%                            particles, and the locations of the particles.
%
% Note: For best performance the sum of all values in the BeadFilter should
% be zero. This can be achieved by removing the mean BeadFilter value from 
% each value. This ensures that any background intensities are ignored.

%% Locate beads in a still image frame
%
% The coordinate system used is a MATLAB matrix coordinate system with %X=Y=0%
% in the top left corner of the image! (Row, Column)


%% General settings
ROI_Size = size(BeadImage);


%% Coarse location of beads
% Construct an image of a point source of light with no pixel offset



% Convolute with image
ImageFrameConvolution = conv2(double(ImageFrame), double(BeadImage), 'same');
% Find local maxima
LocalPeaks = imregionalmax(ImageFrameConvolution, 8);

% Determine the strongest correlation between the constructed bead image
% and the observed frame image. The correlation between the
% constructed bead image and the image has to be EdgeDetectionLimit *
% MaxConvolution to be regarded as a bead
MaxConvolution = max(max(ImageFrameConvolution));
%%%%%%%%%%%%%%%%%%%%%%% MinimumCorrelationFraction = 0.2;
% A similar check is done locally:
%%%%%%%%%%%%%%%%%%%%%%%MinimumLocalCorrelationFraction = 0.5;
% Find beads and count the number of found ones
% FindBeads();
% % Roughly count the number of beads in the frame
% N_Beads = sum(sum(Beads));
% % Store bead locations in matrix
% Coordinates = zeros(N_Beads, 2);


% ROI radius
ROI_R = floor((ROI_Size - 1) / 2);

% Prepare the mask for the local maximum test
[X, Y] = meshgrid(1:ROI_Size(1), 1:ROI_Size(2));
X = X - ROI_R(1) - 1;
Y = Y - ROI_R(2) - 1;
LocalMaximumMask = sqrt(X.^2 + Y.^2) < ROI_R(1);
LocalMaximumMask_Ring = sqrt(X.^2 + Y.^2) < ROI_R(1)+1 & sqrt(X.^2 + Y.^2) > ROI_R(1) - 1;

LocalMaximumRingMap = zeros(size(ImageFrameConvolution));
Location = FindBeads();

% Save the number of found objects
Population.N_Object = size(Location,1);
Population.Location = Location;

% Plot the final result if requested

if Diagnostic
    figure;
    h_Zoom = zoom();
    h_Zoom.ActionPostCallback = @ZoomCallBack;
    h_Pan = pan();
    h_Pan.ActionPostCallback = @ZoomCallBack;
    
    % Plot original image
    h_Image = subplot(4,4,[1,2,5,6]);
    PlotImage();
    
    % Plot convolution
    h_Convolution = subplot(4,4,[9,10,13,14]);
    ResultImage = DetermineFoundBeads(0.7, 0.8);
    PlotResult();
    
    % Convolution histogram
    h_ConvolutionHistogram = subplot(4,4,3);
    HistData = LocalPeaks.*ImageFrameConvolution;
    HistData(HistData > 0);
    histogram(h_ConvolutionHistogram, HistData(:) ./ MaxConvolution, 0:0.01:1);
    h_ConvolutionHistogram.XLim = [0 1];
    h_ConvolutionHistogram.YLim = [0.1 10000];
    h_ConvolutionHistogram.YScale = 'log';
    
    SliderHeights = 0.02;
    TextHeight = 0.04;
    
    MaxInt = double(max(max(ImageFrame)));
    
    h_MinImageIntensity = uicontrol('Style', 'slider', 'Max', MaxInt, 'sliderstep', [1./MaxInt 5./MaxInt], 'value', 0, 'Units' ,'normalized', 'Position', [0.8 0.9 0.15 SliderHeights]);
    h_MinImageIntensityText = uicontrol('Style', 'text', 'string', 'Maximum plotted intensity: 0', 'Units' ,'normalized', 'Position', [0.8 0.925 0.15 TextHeight]);
    h_MinImageIntensity.Callback = @Call_IntensityMin;
    
    h_MaxImageIntensity = uicontrol('Style', 'slider', 'Min', 1, 'Max', MaxInt, 'sliderstep', [1./MaxInt 5./MaxInt], 'value', max(max(ImageFrame)), 'Units' ,'normalized', 'Position', [0.8 0.8 0.15 SliderHeights]);
    h_MaxImageIntensityText = uicontrol('Style', 'text', 'string', ['Minimum plotted intensity: ' num2str(max(max(ImageFrame)))], 'Units' ,'normalized', 'Position', [0.8 0.825 0.15 TextHeight]);
    h_MaxImageIntensity.Callback = @Call_IntensityMax;
    
    h_CorrSlider = uicontrol('Style', 'slider', 'sliderstep', [0.01 0.1], 'value', 0.7, 'Units' ,'normalized', 'Position', [0.8 0.7 0.15 SliderHeights]);
    h_CorrSliderText = uicontrol('Style', 'text', 'string', 'Minimum correlation factor: 0.7', 'Units' ,'normalized', 'Position', [0.8 0.725 0.15 TextHeight]);
    h_CorrSlider.Callback = @Call_CorrValue;
    
    h_LocalCorrSlider = uicontrol('Style', 'slider', 'Min', 0.01, 'sliderstep', [0.01 0.1], 'value', 0.8, 'Units' ,'normalized', 'Position', [0.8 0.6 0.15 SliderHeights]);
    h_LocalCorrSliderText = uicontrol('Style', 'text', 'string', 'Maximum local factor: 0.8', 'Units' ,'normalized', 'Position', [0.8 0.625 0.15 TextHeight]);
    h_LocalCorrSlider.Callback = @Call_LocalCorrValue;
    
  
    
    
end

    function Location =  FindBeads()
        % Determine the location of all bead candidates
        Beads = (ImageFrameConvolution .* LocalPeaks) > MinimumCorrelationFraction * MaxConvolution;
        % The running index for the beads
        i_FoundBead = 1;
        
        Coordinates = [];
        % Loop over the whole image and store the bead locations
        for my = 6:(size(ImageFrameConvolution, 1)-5)
            for nx = 6:(size(ImageFrameConvolution, 2)-5)
                if Beads(my,nx)
                    
                    % Determine the ROI boundaries
                    ROI_X_Min = round(nx - ROI_R(1));
                    ROI_X_Max = round(nx + ROI_R(1));
                    ROI_Y_Min = round(my - ROI_R(2));
                    ROI_Y_Max = round(my + ROI_R(2));
                    
                    % Flag to indicate edge objects (ROI intersects edge)
                    Edge = false;
                    % Raise flag if object is at the edge
                    
                    if ROI_X_Min < 1 + EdgeWidth
                        ROI_X_Min = 1;
                        Edge = true;
                    end
                    if ROI_X_Max > size(ImageFrame, 2) - EdgeWidth
                        ROI_X_Max = size(ImageFrame, 2);
                        Edge = true;
                    end
                    if ROI_Y_Min < 1 + EdgeWidth
                        ROI_Y_Min = 1;
                        Edge = true;
                    end
                    if ROI_Y_Max > size(ImageFrame, 1) - EdgeWidth
                        ROI_Y_Max = size(ImageFrame, 1);
                        Edge = true;
                    end
                    
                    % If a bead is close to the edge, ignore it
                    if ~Edge
                        LocalMaximum_Ring = max(max(LocalMaximumMask_Ring.*ImageFrameConvolution(ROI_Y_Min:ROI_Y_Max, ROI_X_Min:ROI_X_Max)));
                        
                        LocalMaximumRingMap(my,nx) = LocalMaximum_Ring;
                        
                        % Check if the local maximum is >
                        % MinimumLocalCorrelationFraction * LocalMaximum;
                        if LocalMaximum_Ring < ImageFrameConvolution(my,nx) * MaximumLocalCorrelationFraction;
                            Coordinates(i_FoundBead, :) = [nx, my];
                            i_FoundBead = i_FoundBead + 1;
                        end
                        
                    end
                    
                end
            end
        end
        N_Object = i_FoundBead - 1;
        Location = Coordinates(1:N_Object, :);
        
    end

%% Callback functions
    function ZoomCallBack(obj,event_obj)
        if event_obj.Axes == h_Convolution || event_obj.Axes == h_Image
            h_Image.XLim = event_obj.Axes.XLim;
            h_Image.YLim = event_obj.Axes.YLim;
            h_Convolution.XLim = event_obj.Axes.XLim;
            h_Convolution.YLim = event_obj.Axes.YLim;
        end
    end

    % Plot frame
    function PlotImage()
        axes(h_Image);
        imagesc(ImageFrame);
        title('Image', 'FontSize', 16);
        CMap = colormap('gray');
        CMap(64,:) = [1 0 0];
        colormap(CMap);
        axis image;
    end

    % Plot diagnostic result
    function PlotResult()
        axes(h_Convolution);
        
        imagesc(ResultImage);
        hold on;
        axis image;
        %        colormap gray;
        xlabel('X coordinate', 'FontSize', 16);
        ylabel('Y coordinate', 'FontSize', 16);
        title('Convolution', 'FontSize', 16);
        
        h_Convolution.XLim = h_Image.XLim;
        h_Convolution.YLim = h_Image.YLim;
        
    end


    % Determine diagnostic result image
    function ResultImage = DetermineFoundBeads(CorrFraction, LocalMaximumCorr)
        Bd = fspecial('disk', ROI_R(1).*2);
        
        % Determine the result
        Result = (ImageFrameConvolution .* LocalPeaks) > CorrFraction * MaxConvolution & ...
            LocalMaximumRingMap < ImageFrameConvolution .* LocalMaximumCorr;
        
        ResultImage = conv2(double(Result), ...
            double(Bd), ...
            'same') > 0;
        ResultImage = (ImageFrameConvolution - min(min(ImageFrameConvolution))) .* ResultImage;

    end

    % Callback: Correlation value
    function Call_CorrValue(varargin)
        h_CorrSlider.Value = round(h_CorrSlider.Value .* 100) ./ 100;
        MinimumCorrelationFraction = h_CorrSlider.Value;        
        FindBeads();
        ResultImage = DetermineFoundBeads(h_CorrSlider.Value, h_LocalCorrSlider.Value);
        h_CorrSliderText.String = ['Minimum correlation factor: ' num2str(h_CorrSlider.Value)];
        disp(['Minimum correlation factor: ' num2str(h_CorrSlider.Value)]);
        PlotResult();
    end
    
    % Callback: Local correlation value
    function Call_LocalCorrValue(varargin)
        h_LocalCorrSlider.Value = round(h_LocalCorrSlider.Value .* 100) ./ 100;
        
        MaximumLocalCorrelationFraction = h_LocalCorrSlider.Value;        
        FindBeads();
        ResultImage = DetermineFoundBeads(h_CorrSlider.Value, h_LocalCorrSlider.Value);
        h_LocalCorrSliderText.String = ['Maximum local correlation factor: ' num2str(h_LocalCorrSlider.Value)];
        disp(['Maximum local correlation factor: ' num2str(h_LocalCorrSlider.Value)]);
        PlotResult();
    end

    function Call_IntensityMin(varargin)
        h_MinImageIntensity.Value = round(h_MinImageIntensity.Value);
        if h_MinImageIntensity.Value >= h_MaxImageIntensity.Value
            h_MinImageIntensity.Value =  h_MaxImageIntensity.Value - 1;
        end
        h_Image.CLim = [h_MinImageIntensity.Value h_MaxImageIntensity.Value];
%         h_Convolution.CLim = [h_MinImageIntensity.Value h_MaxImageIntensity.Value];
        h_MinImageIntensityText.String = ['Minimum plotted intensity: ' num2str(h_MinImageIntensity.Value)];
%         disp(['Minimum plotted intensity: ' num2str(h_MinImageIntensity.Value)]);
    end

    function Call_IntensityMax(varargin)
        h_MaxImageIntensity.Value = round(h_MaxImageIntensity.Value);
        if h_MinImageIntensity.Value >= h_MaxImageIntensity.Value
            h_MaxImageIntensity.Value =  h_MinImageIntensity.Value + 1;
        end
        h_Image.CLim = [h_MinImageIntensity.Value h_MaxImageIntensity.Value];
%         h_Convolution.CLim = [h_MinImageIntensity.Value h_MaxImageIntensity.Value];
        h_MaxImageIntensityText.String = ['Maximum plotted intensity: ' num2str(h_MaxImageIntensity.Value)];
%         disp(['Maximum plotted intensity: ' num2str(h_MaxImageIntensity.Value)]);
    end

end


