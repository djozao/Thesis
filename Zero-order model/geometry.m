function object = geometry(CR, L, D, A, Emi)

    if A>D
        fprintf('Aperture diameter cannot be bigger than the diameter itself')
    end
    if(CR == 0)
        object = cylinder(D,A,L,Emi);
    elseif CR == 1
        object = cone(D,A,L, Emi);
    end
end