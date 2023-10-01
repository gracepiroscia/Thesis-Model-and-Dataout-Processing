function sem_out = SEM(Array)
%SEM Calculate the standard error on the mean

    std_array = std(Array);
    n = length(Array);
    sem_out = std_array/sqrt(n);
end

