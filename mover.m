function [nuevoA, nuevoB,vida_new,duplicar_new,contador_new] = mover(A,B,m_vida,m_duplicar,m_contador)
        n = length(A);
        % Buscamos las posiciones de las células CAR-T
        [filaM,columnaM] = find(A==1);
        % Recorremos todas las CAR-T para moverlas
        for j = 1:length(filaM)
            % Si la posición correspondiente en m_duplicar es 0 quiere decir que
            % se puede mover, si no está igual a 0 es que la CAR-T está en
            % fase de duplicacion  y no se mueve.
            if m_duplicar(filaM(j),columnaM(j)) == 0
                % Para inicializar los vectores donde guardamos las posiciones
                % a las que se puede mover una CAR-T.
                F = zeros(1,2);
                C = zeros(1,2);
                % Definimos contador para calcular las probabilidades de moverse a cada posición.
                count = 0;
                % Cuadrado interior de la malla
                if filaM(j) > 1 & columnaM(j) > 1 & filaM(j) < n & columnaM(j) < n
                    if A(filaM(j)-1,columnaM(j)-1) == 0 % Comprobamos que la casilla está vacía, si no está vacía no se puede mover ahí.
                        count = count + 1; % Sumamos uno al contador para calcular luego las probabilidades
                        F(count) = filaM(j)-1; % Guardamos la fila
                        C(count) = columnaM(j) -1; % Guardamos la columna
                    end
                    if A(filaM(j)-1,columnaM(j)) == 0
                        count = count + 1;
                        F(count) = filaM(j)-1;
                        C(count) = columnaM(j);
                    end
                    if A(filaM(j)-1,columnaM(j)+1) == 0
                        count = count + 1;
                        F(count) = filaM(j)-1;
                        C(count) = columnaM(j) +1;
                    end
                    if A(filaM(j),columnaM(j)-1) == 0
                        count = count + 1;
                        F(count) = filaM(j);
                        C(count) = columnaM(j) -1;
                    end
                    if A(filaM(j),columnaM(j)+1) == 0
                        count = count + 1;
                        F(count) = filaM(j);
                        C(count) = columnaM(j) +1;
                    end
                    if A(filaM(j)+1,columnaM(j)-1) == 0
                        count = count + 1;
                        F(count) = filaM(j)+1;
                        C(count) = columnaM(j) -1;
                    end
                    if A(filaM(j)+1,columnaM(j)) == 0
                        count = count + 1;
                        F(count) = filaM(j)+1;
                        C(count) = columnaM(j) ;
                    end
                    if A(filaM(j)+1,columnaM(j)+1) == 0
                        count = count + 1;
                        F(count) = filaM(j)+1;
                        C(count) = columnaM(j) + 1;
                    end
                    if count > 0 % Si el contador es >0 es que se puede mover a algún sitio, sino se queda donde está.
                        k = round((count-1)*rand(1,1))+1; % Sacamos un num aleatorio entre 1 y count y esa es la posición a la que se mueve.
                        % la posición a la que se mueve la CAR-T es (F(k),C(k))
                        A(filaM(j), columnaM(j)) = 0; % Ponemos a 0 la posición inicial de la CAR-T
                        A(F(k),C(k)) = 1; % Ponemos a 1 la nueva posición de la CAR-T
                        % Actualizamos la matriz vida
                        m_vida(F(k),C(k)) = m_vida(filaM(j), columnaM(j));
                        m_vida(filaM(j), columnaM(j)) = 0;
                        % Actualizamos la matriz contador
                        m_contador(F(k),C(k)) = m_contador(filaM(j), columnaM(j));
                        m_contador(filaM(j), columnaM(j)) = 0;
                        % Tenemos que cambiar también las posiciones de las CAR-T en la matriz B, en esta matriz se guarda la
                        % carga que tiene cada CAR-T disponible.
                        B(F(k),C(k)) = B(filaM(j), columnaM(j));
                        B(filaM(j), columnaM(j)) = 0;
                    end
                end
                % Esquina superior izquierda de la malla
                count = 0;
                if filaM(j) == 1 & columnaM(j) == 1
                    if A(filaM(j),columnaM(j)+1) == 0
                        count = count + 1;
                        F(count) = filaM(j);
                        C(count) = columnaM(j) +1;
                    end
                    if A(filaM(j)+1,columnaM(j)) == 0
                        count = count + 1;
                        F(count) = filaM(j)+1;
                        C(count) = columnaM(j);
                    end
                    if A(filaM(j)+1,columnaM(j)+1) == 0
                        count = count + 1;
                        F(count) = filaM(j)+1;
                        C(count) = columnaM(j) +1;
                    end
                    if count > 0
                        k = round((count-1)*rand(1,1))+1;
                        A(filaM(j), columnaM(j)) = 0;
                        A(F(k),C(k)) = 1;
                        m_vida(F(k),C(k)) = m_vida(filaM(j), columnaM(j));
                        m_vida(filaM(j), columnaM(j)) = 0;
                        m_contador(F(k),C(k)) = m_contador(filaM(j), columnaM(j));
                        m_contador(filaM(j), columnaM(j)) = 0;
                        B(F(k),C(k)) = B(filaM(j), columnaM(j));
                        B(filaM(j), columnaM(j)) = 0;
                    end
                end
                % Esquina superior derecha de la malla
                count = 0;
                if filaM(j) == 1 & columnaM(j) == n
                    if A(filaM(j),columnaM(j)-1) == 0
                        count = count + 1;
                        F(count) = filaM(j);
                        C(count) = columnaM(j) - 1;
                    end
                    if A(filaM(j)+1,columnaM(j)-1) == 0
                        count = count + 1;
                        F(count) = filaM(j)+1;
                        C(count) = columnaM(j)-1;
                    end
                    if A(filaM(j)+1,columnaM(j)) == 0
                        count = count + 1;
                        F(count) = filaM(j)+1;
                        C(count) = columnaM(j);
                    end
                    if count > 0
                        k = round((count-1)*rand(1,1))+1;
                        A(filaM(j), columnaM(j)) = 0;
                        A(F(k),C(k)) = 1;
                        m_vida(F(k),C(k)) = m_vida(filaM(j), columnaM(j));
                        m_vida(filaM(j), columnaM(j)) = 0;
                        m_contador(F(k),C(k)) = m_contador(filaM(j), columnaM(j));
                        m_contador(filaM(j), columnaM(j)) = 0;
                        B(F(k),C(k)) = B(filaM(j), columnaM(j));
                        B(filaM(j), columnaM(j)) = 0;
                    end
                end
                % Esquina inferior izquierda de la malla
                count = 0;
                if filaM(j) == n & columnaM(j) == 1
                    if A(filaM(j)-1,columnaM(j)) == 0
                        count = count + 1;
                        F(count) = filaM(j)-1;
                        C(count) = columnaM(j);
                    end
                    if A(filaM(j)-1,columnaM(j)+1) == 0
                        count = count + 1;
                        F(count) = filaM(j)-1;
                        C(count) = columnaM(j)+1;
                    end
                    if A(filaM(j),columnaM(j)+1) == 0
                        count = count + 1;
                        F(count) = filaM(j);
                        C(count) = columnaM(j)+1;
                    end
                    if count > 0
                        k = round((count-1)*rand(1,1))+1;
                        A(filaM(j), columnaM(j)) = 0;
                        A(F(k),C(k)) = 1;
                        m_vida(F(k),C(k)) = m_vida(filaM(j), columnaM(j));
                        m_vida(filaM(j), columnaM(j)) = 0;
                        m_contador(F(k),C(k)) = m_contador(filaM(j), columnaM(j));
                        m_contador(filaM(j), columnaM(j)) = 0;
                        B(F(k),C(k)) = B(filaM(j), columnaM(j));
                        B(filaM(j), columnaM(j)) = 0;
                    end
                end
                % Esquina inferior derecha de la malla
                count = 0;
                if filaM(j) == n & columnaM(j) == n
                    if A(filaM(j)-1,columnaM(j)-1) == 0
                        count = count + 1;
                        F(count) = filaM(j)-1;
                        C(count) = columnaM(j)-1;
                    end
                    if A(filaM(j)-1,columnaM(j)) == 0
                        count = count + 1;
                        F(count) = filaM(j)-1;
                        C(count) = columnaM(j);
                    end
                    if A(filaM(j),columnaM(j)-1) == 0
                        count = count + 1;
                        F(count) = filaM(j);
                        C(count) = columnaM(j)-1;
                    end
                    if count > 0
                        k = round((count-1)*rand(1,1))+1;
                        A(filaM(j), columnaM(j)) = 0;
                        A(F(k),C(k)) = 1;
                        m_vida(F(k),C(k)) = m_vida(filaM(j), columnaM(j));
                        m_vida(filaM(j), columnaM(j)) = 0;
                        m_contador(F(k),C(k)) = m_contador(filaM(j), columnaM(j));
                        m_contador(filaM(j), columnaM(j)) = 0;
                        B(F(k),C(k)) = B(filaM(j), columnaM(j));
                        B(filaM(j), columnaM(j)) = 0;
                    end
                end
                % Borde superior de la malla
                count = 0;
                if filaM(j) == 1 & columnaM(j) > 1 & columnaM(j) < n
                    if A(filaM(j),columnaM(j)-1) == 0
                        count = count + 1;
                        F(count) = filaM(j);
                        C(count) = columnaM(j) -1;
                    end
                    if A(filaM(j),columnaM(j)+1) == 0
                        count = count + 1;
                        F(count) = filaM(j);
                        C(count) = columnaM(j) +1;
                    end
                    if A(filaM(j)+1,columnaM(j)-1) == 0
                        count = count + 1;
                        F(count) = filaM(j)+1;
                        C(count) = columnaM(j) -1;
                    end
                    if A(filaM(j)+1,columnaM(j)) == 0
                        count = count + 1;
                        F(count) = filaM(j)+1;
                        C(count) = columnaM(j) ;
                    end
                    if A(filaM(j)+1,columnaM(j)+1) == 0
                        count = count + 1;
                        F(count) = filaM(j)+1;
                        C(count) = columnaM(j) + 1;
                    end
                    if count > 0
                        k = round((count-1)*rand(1,1))+1;
                        A(filaM(j), columnaM(j)) = 0;
                        A(F(k),C(k)) = 1;
                        m_vida(F(k),C(k)) = m_vida(filaM(j), columnaM(j));
                        m_vida(filaM(j), columnaM(j)) = 0;
                        m_contador(F(k),C(k)) = m_contador(filaM(j), columnaM(j));
                        m_contador(filaM(j), columnaM(j)) = 0;
                        B(F(k),C(k)) = B(filaM(j), columnaM(j));
                        B(filaM(j), columnaM(j)) = 0;
                    end
                end
                % Borde inferior de la malla
                count = 0;
                if filaM(j) == n & columnaM(j) > 1 & columnaM(j) < n
                    if A(filaM(j)-1,columnaM(j)-1) == 0
                        count = count + 1;
                        F(count) = filaM(j)-1;
                        C(count) = columnaM(j) -1;
                    end
                    if A(filaM(j)-1,columnaM(j)) == 0
                        count = count + 1;
                        F(count) = filaM(j)-1;
                        C(count) = columnaM(j);
                    end
                    if A(filaM(j)-1,columnaM(j)+1) == 0
                        count = count + 1;
                        F(count) = filaM(j)-1;
                        C(count) = columnaM(j) +1;
                    end
                    if A(filaM(j),columnaM(j)-1) == 0
                        count = count + 1;
                        F(count) = filaM(j);
                        C(count) = columnaM(j) -1;
                    end
                    if A(filaM(j),columnaM(j)+1) == 0
                        count = count + 1;
                        F(count) = filaM(j);
                        C(count) = columnaM(j) +1;
                    end
                    if count > 0
                        k = round((count-1)*rand(1,1))+1;
                        A(filaM(j), columnaM(j)) = 0;
                        A(F(k),C(k)) = 1;
                        m_vida(F(k),C(k)) = m_vida(filaM(j), columnaM(j));
                        m_vida(filaM(j), columnaM(j)) = 0;
                        m_contador(F(k),C(k)) = m_contador(filaM(j), columnaM(j));
                        m_contador(filaM(j), columnaM(j)) = 0;
                        B(F(k),C(k)) = B(filaM(j), columnaM(j));
                        B(filaM(j), columnaM(j)) = 0;
                    end
                end
                % Borde izquierdo de la malla
                count = 0;
                if columnaM(j) == 1 & filaM(j) < n & filaM(j) > 1
                    if A(filaM(j)-1,columnaM(j)) == 0
                        count = count + 1;
                        F(count) = filaM(j)-1;
                        C(count) = columnaM(j);
                    end
                    if A(filaM(j)-1,columnaM(j)+1) == 0
                        count = count + 1;
                        F(count) = filaM(j)-1;
                        C(count) = columnaM(j) +1;
                    end
                    if A(filaM(j),columnaM(j)+1) == 0
                        count = count + 1;
                        F(count) = filaM(j);
                        C(count) = columnaM(j) +1;
                    end
                    if A(filaM(j)+1,columnaM(j)) == 0
                        count = count + 1;
                        F(count) = filaM(j)+1;
                        C(count) = columnaM(j) ;
                    end
                    if A(filaM(j)+1,columnaM(j)+1) == 0
                        count = count + 1;
                        F(count) = filaM(j) + 1;
                        C(count) = columnaM(j) + 1;
                    end
                    if count > 0
                        k = round((count-1)*rand(1,1))+1;
                        A(filaM(j), columnaM(j)) = 0;
                        A(F(k),C(k)) = 1;
                        m_vida(F(k),C(k)) = m_vida(filaM(j), columnaM(j));
                        m_vida(filaM(j), columnaM(j)) = 0;
                        m_contador(F(k),C(k)) = m_contador(filaM(j), columnaM(j));
                        m_contador(filaM(j), columnaM(j)) = 0;
                        B(F(k),C(k)) = B(filaM(j), columnaM(j));
                        B(filaM(j), columnaM(j)) = 0;
                    end
                end
                % Borde derecho de la malla
                count = 0;
                if columnaM(j) == n & filaM(j) < n & filaM(j) > 1
                    if A(filaM(j)-1,columnaM(j)-1) == 0
                        count = count + 1;
                        F(count) = filaM(j)-1;
                        C(count) = columnaM(j) -1;
                    end
                    if A(filaM(j)-1,columnaM(j)) == 0
                        count = count + 1;
                        F(count) = filaM(j)-1;
                        C(count) = columnaM(j);
                    end
                    if A(filaM(j),columnaM(j)-1) == 0
                        count = count + 1;
                        F(count) = filaM(j);
                        C(count) = columnaM(j) -1;
                    end
                    if A(filaM(j)+1,columnaM(j)-1) == 0
                        count = count + 1;
                        F(count) = filaM(j)+1;
                        C(count) = columnaM(j) -1;
                    end
                    if A(filaM(j)+1,columnaM(j)) == 0
                        count = count + 1;
                        F(count) = filaM(j)+1;
                        C(count) = columnaM(j) ;
                    end
                    if count > 0
                        k = round((count-1)*rand(1,1))+1;
                        A(filaM(j), columnaM(j)) = 0;
                        A(F(k),C(k)) = 1;
                        m_vida(F(k),C(k)) = m_vida(filaM(j), columnaM(j));
                        m_vida(filaM(j), columnaM(j)) = 0;
                        m_contador(F(k),C(k)) = m_contador(filaM(j), columnaM(j));
                        m_contador(filaM(j), columnaM(j)) = 0;
                        B(F(k),C(k)) = B(filaM(j), columnaM(j));
                        B(filaM(j), columnaM(j)) = 0;
                    end
                end
            end
        end
    nuevoA = A;
    nuevoB = B;
    vida_new = m_vida;
    duplicar_new = m_duplicar;
    contador_new = m_contador;
end