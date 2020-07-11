function [HiPass, LoPass] = BackgroundCorrection_WaveletSet(FrameData, Level)
% Function BackgroundCorrection_WaveletSet
%   Image Analysis Framework 
%   Version 1.2 (Emiel Visser)
% This function splits the background and the signal in an image based on
% wavelet filtering of the image.
% 
% *Input parameters:
%         FrameData ::    The image to be filtered (MxN matrix)
%         Level ::        The wavelet level seperating between the low and
%                           high pass
%     Note: A _higher_ 'Level' parameter will _lower_ the frequency limit
%     in the image that is used to seperate the background and signal
%     Thus, a lower level will reject more background. See what works for
%     your data.
%
%
% *Output:
%         HiPass ::       Image component above the wavelet level
%         LoPass ::       Image component below the wavelet level
% 

% We do a correction using a wavelet filter employing a set of wavelets
% Hard-coded selection of wavelets


Wavelets = {'bior3.3','bior3.5','bior3.7','bior3.9','bior4.4','bior5.5','bior6.8'};

% Prepare lo-pass data
LoPass = zeros(size(FrameData));

% Cycle through every wavelet, add contribution to background
for i_ndex = 1:numel(Wavelets)
    [C,S] = wavedec2(FrameData, Level, Wavelets{i_ndex});
    LoPass = LoPass + wrcoef2('a',C,S,Wavelets{i_ndex}, Level);
end

% Normalize due to number of wavelets
LoPass = LoPass ./ numel(Wavelets);
% Determine hi-pass by subtracting lo-pass component
HiPass = double(FrameData) - LoPass;
end