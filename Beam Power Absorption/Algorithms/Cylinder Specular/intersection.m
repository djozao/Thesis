function intersection = intersection(A,B)
    if B(1) == B(2)
        intersection = 0;
        return;
    end
    if A(1)> A(2)
        aux = A(2);
        A(2) = A(1);
        A(1) = aux;
    end
    if B(1)> B(2)
        aux = B(2);
        B(2) = B(1);
        B(1) = aux;
    end

    if min(A(1),B(1)) == A(1)
        if min(B(1),A(2)) == A(2)
            intersection = 0;
            return;
        elseif min(B(1), A(2)) == B(1)
            if min(B(2), A(2)) == B(2)
                intersection = [B(1), B(2)];
                return;
            else
                intersection = [B(1), A(2)];
                return;
            end
        end
    elseif min(A(1), B(1)) == B(1)
        if min(B(2), A(1)) == B(2)
            intersection = 0;
            return;
        elseif min(B(2), A(1)) == A(1)
            if min(A(2), B(2)) == B(2)
                intersection = [A(1), B(2)];
                return;
            else 
                intersection = [A(1), A(2)];
                return
            end
        end
    end       

end