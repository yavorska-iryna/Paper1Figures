function SEM = sem(data)
SEM = nanstd(data)/sqrt(length(data));
end
