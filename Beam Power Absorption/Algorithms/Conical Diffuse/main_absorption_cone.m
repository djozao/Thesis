Results = zeros(5*10*10*17,5);
n1=1; % surface top nodes
n2=1; % lateral surface nodes
i = 1;
D = 1;
for alpha = [5, 15, 25, 35, 45]
    for abs = [0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8,0.9,1]
        for AD = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 ,0.9, 1]
            for LD = [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 3, 4, 5, 6, 7, 8, 9, 10]
                Results(i,1) = alpha;
                Results(i,2) = abs;
                Results(i,3) = AD;
                Results(i,4) = LD;
                Results(i,5) = cone_fem_absorbed(alpha,LD,D,AD,abs,n1,n2);
                i=i+1
            end
        end
    end
end