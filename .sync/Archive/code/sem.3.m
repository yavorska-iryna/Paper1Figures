function SEM = sem(data)
SEM = std(data)/sqrt(size(data));
end
