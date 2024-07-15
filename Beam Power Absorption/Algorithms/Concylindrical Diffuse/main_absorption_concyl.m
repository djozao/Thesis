Results = zeros(5*10*10*17*4,6);
n1=100; % surface top nodes
n2=100; % lateral surface nodes
n3=100; % bottom surface nodes
i = 1;
D = 1;
for alpha = [5, 25, 45]
    for abs = [0.1, 0.3, 0.5, 0.7, 0.9]
        for AD = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 ,0.9, 1]
            for LD = [0.5, 0.75, 1, 1.5, 2, 2.5, 4, 7, 10]
                for CR = [0.1, 0.25, 0.5, 0.75]
                Results(i,1) = alpha;
                Results(i,2) = abs;
                Results(i,3) = AD;
                Results(i,4) = LD;
                Results(i,5) = CR;
                Results(i,6) =  concyl_fem_absorbed(alpha,LD,D,AD, CR, abs,n1,n2,n3);
                i=i+1
            
                end
            end
        end
    end
end