function Gebhart_matrix = Gebhart_3D_Cylinder(L, D, A, emi, n_1, n_2, n_3)


vf_matrix = zeros(1+n_1+n_2+n_3, 1+n_1+n_2+n_3);
vf_discdisc = @(R2,R1,h) 1./2.*(1+ (1+(R2/h)^2)/((R1/h)^2)-sqrt((1+ (1+(R2/h)^2)/((R1/h)^2)).^2-4*(R2/R1)^2));
vf_stripstrip = @(R,h) 1 - 1/2*(sqrt(4*((R/h)^2)+1)-1)/((R/h));

for i = 1:n_3
    %F_3i_0
    if i == 1
        vf_matrix(1+n_1+n_2+i, 1) = vf_discdisc(A/2, D/n_3/2, L);
    else
        vf_matrix(1+n_1+n_2+i, 1) = ((D/n_3/2*i)^2*vf_discdisc(A/2, D/n_3/2*i, L)-(D/n_3/2*(i-1))^2*vf_discdisc(A/2, D/n_3/2*(i-1), L))/((D/n_3/2*i)^2-(D/n_3/2*(i-1))^2);
    end
        vf_matrix(1, 1+n_1+n_2+i) = vf_matrix(1+n_1+n_2+i, 1)*((D/n_3/2*i)^2-(D/n_3/2*(i-1))^2)/((A/2)^2);
    %F_3i_1j
    for j = 1:n_1
        if i ==1
            vf_matrix(1+n_1+n_2+i, 1+j) = vf_discdisc(A/2+(D-A)/2/n_1*j, D/n_3/2, L)-vf_discdisc(A/2+(D-A)/2/n_1*(j-1), D/n_3/2, L);
        else
            vf_matrix(1+n_1+n_2+i, 1+j) = (((i*D/n_3/2)^2)*(vf_discdisc(A/2+(D-A)/2/n_1*j, i*D/n_3/2, L)-vf_discdisc(A/2+(D-A)/2/n_1*(j-1), i*D/n_3/2, L)) - ((D/n_3/2*(i-1))^2)*(vf_discdisc(A/2+(D-A)/2/n_1*j, (i-1)*D/n_3/2, L)-vf_discdisc(A/2+(D-A)/2/n_1*(j-1), (i-1)*D/n_3/2, L)))/(((i*D/n_3/2)^2)-(D/n_3/2*(i-1))^2);
        end
        if A == 1
            vf_matrix(1+j, 1+n_1+n_2+i) = 0;
        else
            vf_matrix(1+j, 1+n_1+n_2+i) = vf_matrix(1+n_1+n_2+i, 1+j)*((D/n_3/2*i)^2-(D/n_3/2*(i-1))^2)/((A/2+(D-A)/2/n_1*j)^2-(A/2+(D-A)/2/n_1*(j-1))^2);        
        end
    end

    %F_3i_2j
    for j = 1:n_2
        if i == 1
            if j  ==1
                vf_matrix(1+n_1+n_2+i, 1+n_1+j) = 1-vf_discdisc(D/2, D/n_3/2, L/n_2*j);
            else
                vf_matrix(1+n_1+n_2+i, 1+n_1+j) = vf_discdisc(D/2, D/n_3/2, L/n_2*(j-1))-vf_discdisc(D/2, D/n_3/2, L/n_2*j);
            end
        else
            if j ==1
                vf_matrix(1+n_1+n_2+i, 1+n_1+j) = (((D/n_3/2*i)^2)*(1-vf_discdisc(D/2, D/n_3/2*i, L/n_2*j))-((D/n_3/2*(i-1))^2)*(1-vf_discdisc(D/2, D/n_3/2*(i-1), L/n_2*j)))/(((i*D/n_3/2)^2)-(D/n_3/2*(i-1))^2);
            else
                vf_matrix(1+n_1+n_2+i, 1+n_1+j) = (((D/n_3/2*i)^2)*(vf_discdisc(D/2, D/n_3/2*i, L/n_2*(j-1))-vf_discdisc(D/2, D/n_3/2*i, L/n_2*j))-((D/n_3/2*(i-1))^2)*(vf_discdisc(D/2, D/n_3/2*(i-1), L/n_2*(j-1))-vf_discdisc(D/2, D/n_3/2*(i-1), L/n_2*j)))/(((i*D/n_3/2)^2)-(D/n_3/2*(i-1))^2);
            end
        end
            vf_matrix(1+n_1+j, 1+n_1+n_2+i) = vf_matrix(1+n_1+n_2+i, 1+n_1+j)*((D/n_3/2*i)^2-(D/n_3/2*(i-1))^2)/(L/n_2*D);
    end
end

% F2

for i = 1:n_2
    %F2i_2i
    for j = i:n_2
        if j == i
            vf_matrix(1+n_1+i, 1+n_1+j) = vf_stripstrip(D/2,L/n_2);
        elseif j == i+1
            vf_matrix(1+n_1+i, 1+n_1+j) = D/(4*L/n_2)*(1-vf_discdisc(D/2,D/2,L/n_2*(j-i))-vf_discdisc(D/2,D/2,L/n_2*(j-i))+vf_discdisc(D/2,D/2,L/n_2*(j-i+1)));
        else
            vf_matrix(1+n_1+i, 1+n_1+j) = D/(4*L/n_2)*(vf_discdisc(D/2,D/2, L/n_2*(j-i-1))-vf_discdisc(D/2,D/2,L/n_2*(j-i))-vf_discdisc(D/2,D/2,L/n_2*(j-i))+vf_discdisc(D/2,D/2,L/n_2*(j-i+1)));
        end
            vf_matrix(1+n_1+j, 1+n_1+i) = vf_matrix(1+n_1+i, 1+n_1+j);
    end
    %F2i_0
    if i == n_2
        vf_matrix(1, 1+n_1+i) = 1-vf_discdisc(D/2, A/2, L-L/n_2*(i-1));
    else
        vf_matrix(1, 1+n_1+i) = vf_discdisc(D/2, A/2, L-L/n_2*(i))-vf_discdisc(D/2, A/2, L-L/n_2*(i-1));
    end
        vf_matrix(1+n_1+i, 1) = vf_matrix(1, 1+n_1+i)*((A/2)^2)/(L/n_2*D);
end

% F_1i_2j
for i = 1:n_1
    if A ==1
        for j = 1:n_2
            vf_matrix(1+i,1+n_1+j) = 0;
            vf_matrix(1+n_1+j,1+i) = 0;
        end
    else
        for j = 1:n_2
            if j == n_2
                vf_matrix(1+i,1+n_1+j) = ( (((D-A)/n_1/2*i+A/2)^2)*(1-vf_discdisc(D/2,(D-A)/n_1/2*i+A/2,L-L/n_2*(j-1))) -((D-A)/n_1/2*(i-1)+A/2)^2*((1-vf_discdisc(D/2,(D-A)/n_1/2*(i-1)+A/2,L-L/n_2*(j-1))))     )/(((D-A)/n_1/2*i+A/2)^2-((D-A)/n_1/2*(i-1)+A/2)^2);
            else
                vf_matrix(1+i,1+n_1+j) = ( (((D-A)/n_1/2*i+A/2)^2)*(vf_discdisc(D/2,(D-A)/n_1/2*i+A/2,L-L/n_2*j)-vf_discdisc(D/2,(D-A)/n_1/2*i+A/2,L-L/n_2*(j-1))) -((D-A)/n_1/2*(i-1)+A/2)^2*((vf_discdisc(D/2,(D-A)/n_1/2*(i-1)+A/2,L-L/n_2*j)-vf_discdisc(D/2,(D-A)/n_1/2*(i-1)+A/2,L-L/n_2*(j-1))))     )/(((D-A)/n_1/2*i+A/2)^2-((D-A)/n_1/2*(i-1)+A/2)^2);
            end
                vf_matrix(1+n_1+j,1+i) = vf_matrix(1+i,1+n_1+j)*(((D-A)/n_1/2*i+A/2)^2-((D-A)/n_1/2*(i-1)+A/2)^2)/(L/n_2*D);
        end
    end
end


A_matrix = vf_matrix*(1-emi);
A_matrix(:,1) = A_matrix(:,1)*0;
A_matrix = A_matrix - eye(1+n_1+n_2+n_3);

B = -vf_matrix;
B(:,2:1+n_1+n_2+n_3) = B(:,2:1+n_1+n_2+n_3)*emi;
Gebhart_matrix = linsolve(A_matrix, B);
end 