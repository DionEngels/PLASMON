%% load in data
load v8_Yuyang/Gaussian

bad_sigma = Localizations(Localizations(:,6) < 0 | Localizations(:,7) < 0 | Localizations(:,6) > 5 | Localizations(:,7) > 5,:);

bad_int = Localizations(Localizations(:,5) < 500 | Localizations(:,5) > 30000, :);

max(Localizations(:,9))