function setdif = setdif(A, intersection)
    if intersection == 0
        setdif = A;
        return;
    end
    if intersection(1)< A(1) || intersection (2)> A(2)
        disp('the intersection is not inside the element')
        return;
    end
    if A(1) == intersection(1)
        if intersection(2) < A(2)
            setdif = [intersection(2), A(2)];
            return;
        else
            setdif = 0;
            return;
        end
    else 
        setdif(1,:) = [A(1), intersection(1)];
        if intersection(2) < A(2)
            setdif(2,:) = [intersection(2), A(2)];
            return;
        else
            return;
        end
    end            
end