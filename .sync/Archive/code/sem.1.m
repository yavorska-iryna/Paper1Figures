function SEM = sem(data)
SEM = std(data,2)./sqrt(size(data,2));
end
