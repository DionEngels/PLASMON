%% load in data
load v8_Yuyang/8_1_Gaussian_Loose

bad_sigma = Localizations(Localizations(:,6) < 0 | Localizations(:,7) < 0 | Localizations(:,6) > 4 | Localizations(:,7) > 4,:);

bad_int = Localizations(Localizations(:,5) < 0 | Localizations(:,5) > 2^16-1, :);

max(Localizations(:,9))