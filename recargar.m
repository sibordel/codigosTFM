function [nuevoB] = recargar(A,B,m_duplicar)
% Vamos a considerar que las matrices entren de las mismas dimensiones que a mover
    n = length(A);
    % Buscamos las posiciones de las CAR-T
    [filaT,columnaT] = find(A==1);
    % Recorremos todas las CAR-T para recargarlas
    for j = 1:length(filaT)
        % Si la posición correspondiente en m_duplicar es 0 quiere decir que
        % se puede recargar, si no está igual a 0 es que la CAR-T está en
        % fase de duplicacion y no se recarga. Saldrá recargada cuando se duplique.
        if m_duplicar(filaT(j),columnaT(j)) == 0
            if B(filaT(j), columnaT(j)) < 9
                B(filaT(j), columnaT(j)) = B(filaT(j), columnaT(j))+1;
            end
        end
    end
    nuevoB = B;
end