function Result = reflection_theory(alpha,emi,AD,LD)

D = 1; % diameter 
L = LD*D;
A = AD*D;
origin = A/2/tan(alpha*pi/180);
Boundary = [atan(A/2/(origin+2*L))*180/pi, alpha];
Reflection_matrix = zeros(1000,3);
Reflection_matrix(1,:) = [0, Boundary(1), 1];
i = 1;
l_coef = 1;

while ~isempty(Boundary)
    d_coef = 1;
    new_values = [atan((d_coef*D-A/2)/(origin+2*l_coef*L))*180/pi,atan((d_coef*D+A/2)/(origin+2*l_coef*L))*180/pi];
    while new_values(1)<= alpha && ~isempty(Boundary)
        boundary_size = size(Boundary,1);
        if boundary_size == 1 && Boundary(1,2) - Boundary(1,1)<10^-3            
            Boundary = [];
        end
        new_values = [atan((d_coef*D-A/2)/(origin+2*l_coef*L))*180/pi,atan((d_coef*D+A/2)/(origin+2*l_coef*L))*180/pi];
        j = 1;
        while j <= boundary_size && ~isempty(Boundary)       
            intersect = intersection(Boundary(j,:), new_values);
            if intersect ~= 0
                Reflection_matrix(i+1,:)= [intersect(1), intersect(2),2*l_coef+d_coef-1];
                i = i+1;
                a = setdif(Boundary(j,:), intersect);
                if a == 0                    
                    Boundary(j,:) = [];
                    j= j-1;
                    boundary_size = boundary_size-1;
                elseif size(a,1) == 1
                    Boundary(j,:) = a(1,:);
                elseif size(a,1) == 2
                    Boundary(j,:) = a(1,:);
                    Boundary = [Boundary; a(2,:)];                    
                end

            end
            j= j+1;
        end
        d_coef = d_coef +1;
    end
    l_coef = l_coef +1;
end

    %Result = 0;
    Result2 = 0;
    for j = 1:i
        %Result = Result+(Reflection_matrix(j,2)^2-Reflection_matrix(j,1)^2)*(1-(1-emi)^Reflection_matrix(j,3));
        Result2 = Result2 + ((1-cos(Reflection_matrix(j,2)*pi/180))-(1-cos(Reflection_matrix(j,1)*pi/180)))*(1-(1-emi)^Reflection_matrix(j,3));
    end
    %Result = Result/(alpha^2);
    Result = Result2 /(1-cos(alpha*pi/180));

end









