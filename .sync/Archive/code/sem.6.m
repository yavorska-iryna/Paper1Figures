function SEM = sem(data)
data(isnan(data)) = [];
SEM = std(data)/sqrt(length(data));
end
