function SEM = sem(data)
% calculates standard error of the mean, excludes NaN. ira 04.27.2020

data(isnan(data)) = []; % get rid of NaN
SEM = std(data)/sqrt(length(data));
end
