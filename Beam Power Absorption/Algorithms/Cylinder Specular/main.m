Results = zeros(5*10*10*17,5);
i = 1;

for alpha = [5, 25, 45]
    for abs = [0.1, 0.3, 0.5, 0.7, 0.9]
        for AD = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 ,0.9, 1]
            for LD = [0.5, 0.75, 1, 1.5, 2, 2.5, 4, 7, 10]
                Results(i,1) = alpha;
                Results(i,2) = abs;
                Results(i,3) = AD;
                Results(i,4) = LD;
                Results(i,5) = reflection_theory(alpha,abs,AD,LD)*100;
                i=i+1
            end
        end
    end
end