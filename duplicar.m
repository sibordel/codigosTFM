function [nuevoA, nuevoB, vida_new, duplicar_new, contador_new] = duplicar(A,B,m_vida,m_duplicar,m_contador)
        n = length(A);
        recarga = 9;
        % Buscamos las posiciones de las células CAR-T que estén en la
        % última hora de la fase de duplicación. Estas son las células que
        % en la siguiente hora se duplicarán.
        [filaD,columnaD] = find(m_duplicar==1);
        % Recorremos todas las CAR-T que vayan a duplicarse
        for j = 1:length(filaD)
            % Si la posición correspondiente en m_duplicar es 0 quiere decir que
            % se puede mover, si no está igual a 0 es que la CAR-T está en
            % fase de duplicacion  y no se mueve.
                % Para inicializar los vectores donde guardamos las posiciones
                % a las que se puede duplicar una células CAR-T
                F = zeros(1,2);
                C = zeros(1,2);
                % Definimos contador para calcular las probabilidades de duplicarse a cada posición.
                count = 0;
                % Cuadrado interior de la malla
                if filaD(j) > 1 & columnaD(j) > 1 & filaD(j) < n & columnaD(j) < n
                    if A(filaD(j) - 1, columnaD(j) - 1) == 0 % Comprobamos que la casilla está vacía, si no está vacía no se puede duplicar ahí.
                        count = count + 1; % Sumamos uno al contador para calcular luego las probabilidades
                        F(count) = filaD(j) - 1; % Guardamos la fila
                        C(count) = columnaD(j) - 1; % Guardamos la columna
                    end
                    if A(filaD(j) - 1, columnaD(j)) == 0
                        count = count + 1;
                        F(count) = filaD(j) - 1;
                        C(count) = columnaD(j);
                    end
                    if A(filaD(j) - 1, columnaD(j) + 1) == 0
                        count = count + 1;
                        F(count) = filaD(j) - 1;
                        C(count) = columnaD(j) + 1;
                    end
                    if A(filaD(j), columnaD(j) - 1) == 0
                        count = count + 1;
                        F(count) = filaD(j);
                        C(count) = columnaD(j) - 1;
                    end
                    if A(filaD(j), columnaD(j) + 1) == 0
                        count = count + 1;
                        F(count) = filaD(j);
                        C(count) = columnaD(j) + 1;
                    end
                    if A(filaD(j) + 1,columnaD(j) - 1) == 0
                        count = count + 1;
                        F(count) = filaD(j) + 1;
                        C(count) = columnaD(j) - 1;
                    end
                    if A(filaD(j) + 1, columnaD(j)) == 0
                        count = count + 1;
                        F(count) = filaD(j) + 1;
                        C(count) = columnaD(j);
                    end
                    if A(filaD(j) + 1,columnaD(j) + 1) == 0
                        count = count + 1;
                        F(count) = filaD(j) + 1;
                        C(count) = columnaD(j) + 1;
                    end
                    if count > 0 % Si el contador es >0 es que se puede duplicar a algún sitio, en otro caso, no se puede duplicar.
                        k = round((count-1)*rand(1,1)) + 1; % Sacamos un num aleatorio entre 1 y count y esa es la posición a la que se duplica.
                        % La posición a la que se mueve la CAR-T es (F(k),C(k))
                        A(F(k),C(k)) = 1; % Ponemos a 1 la posición de la nueva CAR-T
                        % Actualizamos la matriz B. Ambas CAR-T se recargan por completo
                        B(filaD(j), columnaD(j)) = recarga;
                        B(F(k),C(k))= recarga;
                        % Vida completa en la nueva CAR-T
                        m_vida(F(k),C(k)) = 336;
                        % Actualizamos m_duplicar ya que la CAR-T original sale de la fase de duplicación.
                        m_duplicar(filaD(j), columnaD(j)) = 0;
                        % Actualizamos la matriz contador. A ambas se le añade un +1 en sus duplicaciones
                        m_contador(F(k),C(k)) = m_contador(filaD(j), columnaD(j)) + 1;
                        m_contador(filaD(j), columnaD(j)) = m_contador(filaD(j), columnaD(j)) + 1;
                    end
                end
                % Esquina superior izquierda de la malla
                count = 0;
                if filaD(j) == 1 & columnaD(j) == 1
                    if A(filaD(j), columnaD(j) + 1) == 0
                        count = count + 1;
                        F(count) = filaD(j);
                        C(count) = columnaD(j) + 1;
                    end
                    if A(filaD(j) + 1, columnaD(j)) == 0
                        count = count + 1;
                        F(count) = filaD(j) + 1;
                        C(count) = columnaD(j);
                    end
                    if A(filaD(j) + 1, columnaD(j) + 1) == 0
                        count = count + 1;
                        F(count) = filaD(j) + 1;
                        C(count) = columnaD(j) + 1;
                    end
                    if count > 0
                        k = round((count-1)*rand(1,1)) + 1;
                        A(F(k),C(k)) = 1; % Ponemos a 1 la nueva posición de la nueva CAR-T
                        % Actualizamos la matriz B. Ambas CAR-T se recargan
                        B(filaD(j), columnaD(j)) = recarga;
                        B(F(k),C(k))= recarga;
                        % Vida completa en la nueva CAR-T
                        m_vida(F(k),C(k)) = 336;
                        % Actualizamos m_duplicar 
                        m_duplicar(filaD(j), columnaD(j)) = 0;
                        % Actualizamos la matriz contador. A ambas se le añade un +1 en sus duplicaciones
                        m_contador(F(k),C(k)) = m_contador(filaD(j), columnaD(j)) + 1;
                        m_contador(filaD(j), columnaD(j)) = m_contador(filaD(j), columnaD(j)) + 1;
                    end
                end
                % Esquina superior derecha de la malla
                count = 0;
                if filaD(j) == 1 & columnaD(j) == n
                    if A(filaD(j), columnaD(j) - 1) == 0
                        count = count + 1;
                        F(count) = filaD(j);
                        C(count) = columnaD(j) - 1;
                    end
                    if A(filaD(j) + 1, columnaD(j) - 1) == 0
                        count = count + 1;
                        F(count) = filaD(j) + 1;
                        C(count) = columnaD(j) - 1;
                    end
                    if A(filaD(j) + 1, columnaD(j)) == 0
                        count = count + 1;
                        F(count) = filaD(j) + 1;
                        C(count) = columnaD(j);
                    end
                    if count > 0
                        k = round((count-1)*rand(1,1)) + 1;
                        A(F(k),C(k)) = 1; % Ponemos a 1 la nueva posición de la nueva CAR-T
                        % Actualizamos la matriz B. Ambas CAR-T se recargan
                        B(filaD(j), columnaD(j)) = recarga;
                        B(F(k),C(k))= recarga;
                        % Vida completa en la nueva CAR-T
                        m_vida(F(k),C(k)) = 336;
                        % Actualizamos m_duplicar
                        m_duplicar(filaD(j), columnaD(j)) = 0;
                        % Actualizamos la matriz contador. A ambas se le añade un +1 en sus duplicaciones
                        m_contador(F(k),C(k)) = m_contador(filaD(j), columnaD(j)) + 1;
                        m_contador(filaD(j), columnaD(j)) = m_contador(filaD(j), columnaD(j)) + 1;
                    end
                end
                % Esquina inferior izquierda de la malla
                count = 0;
                if filaD(j) == n & columnaD(j) == 1
                    if A(filaD(j) - 1, columnaD(j)) == 0
                        count = count + 1;
                        F(count) = filaD(j) - 1;
                        C(count) = columnaD(j);
                    end
                    if A(filaD(j) - 1, columnaD(j) + 1) == 0
                        count = count + 1;
                        F(count) = filaD(j) - 1;
                        C(count) = columnaD(j) + 1;
                    end
                    if A(filaD(j), columnaD(j) + 1) == 0
                        count = count + 1;
                        F(count) = filaD(j);
                        C(count) = columnaD(j) + 1;
                    end
                    if count > 0
                        k = round((count-1)*rand(1,1)) + 1;
                        A(F(k),C(k)) = 1; % Ponemos a 1 la nueva posición de la nueva CAR-T
                        % Actualizamos la matriz B. Ambas CAR-T se recargan
                        B(filaD(j), columnaD(j)) = recarga;
                        B(F(k),C(k))= recarga;
                        % Vida completa en la nueva CAR-T
                        m_vida(F(k),C(k)) = 336;
                        % Actualizamos m_duplicar
                        m_duplicar(filaD(j), columnaD(j)) = 0;
                        % Actualizamos la matriz contador. A ambas se le añade un +1 en sus duplicaciones
                        m_contador(F(k),C(k)) = m_contador(filaD(j), columnaD(j)) + 1;
                        m_contador(filaD(j), columnaD(j)) = m_contador(filaD(j), columnaD(j)) + 1;
                    end
                end
                % Esquina inferior derecha de la malla
                count = 0;
                if filaD(j) == n & columnaD(j) == n
                    if A(filaD(j) - 1, columnaD(j) - 1) == 0
                        count = count + 1;
                        F(count) = filaD(j) - 1;
                        C(count) = columnaD(j) - 1;
                    end
                    if A(filaD(j) - 1, columnaD(j)) == 0
                        count = count + 1;
                        F(count) = filaD(j) - 1;
                        C(count) = columnaD(j);
                    end
                    if A(filaD(j), columnaD(j) - 1) == 0
                        count = count + 1;
                        F(count) = filaD(j);
                        C(count) = columnaD(j) - 1;
                    end
                    if count > 0
                        k = round((count-1)*rand(1,1)) + 1;
                        A(F(k),C(k)) = 1; % Ponemos a 1 la nueva posición de la nueva CAR-T
                        % Actualizamos la matriz B. Ambas CAR-T se recargan
                        B(filaD(j), columnaD(j)) = recarga;
                        B(F(k),C(k))= recarga;
                        % Vida completa en la nueva CAR-T
                        m_vida(F(k),C(k)) = 336;
                        % Actualizamos m_duplicar
                        m_duplicar(filaD(j), columnaD(j)) = 0;
                        % Actualizamos la matriz contador. A ambas se le añade un +1 en sus duplicaciones
                        m_contador(F(k),C(k)) = m_contador(filaD(j), columnaD(j)) + 1;
                        m_contador(filaD(j), columnaD(j)) = m_contador(filaD(j), columnaD(j)) + 1;
                    end
                end
                % Borde superior de la malla
                count = 0;
                if filaD(j) == 1 & columnaD(j) > 1 & columnaD(j) < n
                    if A(filaD(j), columnaD(j) - 1) == 0
                        count = count + 1;
                        F(count) = filaD(j);
                        C(count) = columnaD(j) - 1;
                    end
                    if A(filaD(j), columnaD(j) + 1) == 0
                        count = count + 1;
                        F(count) = filaD(j);
                        C(count) = columnaD(j) + 1;
                    end
                    if A(filaD(j) + 1, columnaD(j) - 1) == 0
                        count = count + 1;
                        F(count) = filaD(j) + 1;
                        C(count) = columnaD(j) - 1;
                    end
                    if A(filaD(j) + 1, columnaD(j)) == 0
                        count = count + 1;
                        F(count) = filaD(j) + 1;
                        C(count) = columnaD(j);
                    end
                    if A(filaD(j) + 1, columnaD(j) + 1) == 0
                        count = count + 1;
                        F(count) = filaD(j) + 1;
                        C(count) = columnaD(j) + 1;
                    end
                    if count > 0
                        k = round((count-1)*rand(1,1))+1;
                        A(F(k),C(k)) = 1; % Ponemos a 1 la nueva posición de la nueva CAR-T
                        % Actualizamos la matriz B. Ambas CAR-T se recargan
                        B(filaD(j), columnaD(j)) = recarga;
                        B(F(k),C(k))= recarga;
                        % Vida completa en la nueva CAR-T
                        m_vida(F(k),C(k)) = 336;
                        % Actualizamos m_duplicar
                        m_duplicar(filaD(j), columnaD(j)) = 0;
                        % Actualizamos la matriz contador. A ambas se le añade un +1 en sus duplicaciones
                        m_contador(F(k),C(k)) = m_contador(filaD(j), columnaD(j)) + 1;
                        m_contador(filaD(j), columnaD(j)) = m_contador(filaD(j), columnaD(j)) + 1;
                    end
                end
                % Borde inferior de la malla
                count = 0;
                if filaD(j) == n & columnaD(j) > 1 & columnaD(j) < n
                    if A(filaD(j) - 1, columnaD(j) - 1) == 0
                        count = count + 1;
                        F(count) = filaD(j) - 1;
                        C(count) = columnaD(j) - 1;
                    end
                    if A(filaD(j) - 1, columnaD(j)) == 0
                        count = count + 1;
                        F(count) = filaD(j) - 1;
                        C(count) = columnaD(j);
                    end
                    if A(filaD(j) - 1, columnaD(j) + 1) == 0
                        count = count + 1;
                        F(count) = filaD(j) - 1;
                        C(count) = columnaD(j) + 1;
                    end
                    if A(filaD(j), columnaD(j) - 1) == 0
                        count = count + 1;
                        F(count) = filaD(j);
                        C(count) = columnaD(j) - 1;
                    end
                    if A(filaD(j),columnaD(j) + 1) == 0
                        count = count + 1;
                        F(count) = filaD(j);
                        C(count) = columnaD(j) + 1;
                    end
                    if count > 0
                        k = round((count-1)*rand(1,1))+1;
                        A(F(k),C(k)) = 1; % Ponemos a 1 la nueva posición de la nueva CAR-T
                        % Actualizamos la matriz B. Ambas CAR-T se recargan
                        B(filaD(j), columnaD(j)) = recarga;
                        B(F(k),C(k))= recarga;
                        % Vida completa en la nueva CAR-T
                        m_vida(F(k),C(k)) = 336;
                        % Actualizamos m_duplicar
                        m_duplicar(filaD(j), columnaD(j)) = 0;
                        % Actualizamos la matriz contador. A ambas se le añade un +1 en sus duplicaciones
                        m_contador(F(k),C(k)) = m_contador(filaD(j), columnaD(j)) + 1;
                        m_contador(filaD(j), columnaD(j)) = m_contador(filaD(j), columnaD(j)) + 1;
                    end
                end
                % Borde izquierdo de la malla
                count = 0;
                if columnaD(j) == 1 & filaD(j) < n & filaD(j) > 1
                    if A(filaD(j) - 1, columnaD(j)) == 0
                        count = count + 1;
                        F(count) = filaD(j) - 1;
                        C(count) = columnaD(j);
                    end
                    if A(filaD(j) - 1, columnaD(j) + 1) == 0
                        count = count + 1;
                        F(count) = filaD(j) - 1;
                        C(count) = columnaD(j) + 1;
                    end
                    if A(filaD(j), columnaD(j) + 1) == 0
                        count = count + 1;
                        F(count) = filaD(j);
                        C(count) = columnaD(j) + 1;
                    end
                    if A(filaD(j) + 1, columnaD(j)) == 0
                        count = count + 1;
                        F(count) = filaD(j) + 1;
                        C(count) = columnaD(j);
                    end
                    if A(filaD(j) + 1, columnaD(j) + 1) == 0
                        count = count + 1;
                        F(count) = filaD(j) + 1;
                        C(count) = columnaD(j) + 1;
                    end
                    if count > 0
                        k = round((count-1)*rand(1,1))+1;
                        A(F(k),C(k)) = 1; % Ponemos a 1 la nueva posición de la nueva CAR-T
                        % Actualizamos la matriz B. Ambas CAR-T se recargan
                        B(filaD(j), columnaD(j)) = recarga;
                        B(F(k),C(k))= recarga;
                        % Vida completa en la nueva CAR-T
                        m_vida(F(k),C(k)) = 336;
                        % Actualizamos m_duplicar
                        m_duplicar(filaD(j), columnaD(j)) = 0;
                        % Actualizamos la matriz contador. A ambas se le añade un +1 en sus duplicaciones
                        m_contador(F(k),C(k)) = m_contador(filaD(j), columnaD(j)) + 1;
                        m_contador(filaD(j), columnaD(j)) = m_contador(filaD(j), columnaD(j)) + 1;
                    end
                end
                % Borde derecho de la malla
                count = 0;
                if columnaD(j) == n & filaD(j) < n & filaD(j) > 1
                    if A(filaD(j) - 1, columnaD(j) - 1) == 0
                        count = count + 1;
                        F(count) = filaD(j) - 1;
                        C(count) = columnaD(j) - 1;
                    end
                    if A(filaD(j) - 1, columnaD(j)) == 0
                        count = count + 1;
                        F(count) = filaD(j) - 1;
                        C(count) = columnaD(j);
                    end
                    if A(filaD(j), columnaD(j) - 1) == 0
                        count = count + 1;
                        F(count) = filaD(j);
                        C(count) = columnaD(j) - 1;
                    end
                    if A(filaD(j) + 1, columnaD(j) - 1) == 0
                        count = count + 1;
                        F(count) = filaD(j)+1;
                        C(count) = columnaD(j) -1;
                    end
                    if A(filaD(j) + 1, columnaD(j)) == 0
                        count = count + 1;
                        F(count) = filaD(j) + 1;
                        C(count) = columnaD(j);
                    end
                    if count > 0
                        k = round((count-1)*rand(1,1))+1;
                        A(F(k),C(k)) = 1; % Ponemos a 1 la nueva posición de la nueva CAR-T
                        % Actualizamos la matriz B. Ambas CAR-T se recargan
                        B(filaD(j), columnaD(j)) = recarga;
                        B(F(k),C(k))= recarga;
                        % Vida completa en la nueva CAR-T
                        m_vida(F(k),C(k)) = 336;
                        % Actualizamos m_duplicar
                        m_duplicar(filaD(j), columnaD(j)) = 0;
                        % Actualizamos la matriz contador. A ambas se le añade un +1 en sus duplicaciones
                        m_contador(F(k),C(k)) = m_contador(filaD(j), columnaD(j)) + 1;
                        m_contador(filaD(j), columnaD(j)) = m_contador(filaD(j), columnaD(j)) + 1;
                    end
                end
        end
    nuevoA = A;
    nuevoB = B;
    vida_new = m_vida;
    duplicar_new = m_duplicar;
    contador_new = m_contador;
end