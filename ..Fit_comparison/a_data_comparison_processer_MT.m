%clear all
close all
%% load in data
load v7_Yuyang/Last_Fit_7x7
ST = Localizations;
load v7_Yuyang/7_1_Last_Fit_7x7_MT
MT = Localizations;

intensity_ST = ST(ST(:,2) == 1, 5);
intensity_MT = MT(MT(:,2) == 1, 5);