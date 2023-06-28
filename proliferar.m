function [AA,BB,vida_new,duplicar_new,contador_new,FILA,COLUMNA] = proliferar(fila, columna, A, B, m_vida,m_duplicar,m_contador,rho)
    % Meter A y B sin borde en el principal. Entran matrices de tamaño nxn
    n = length(A);
    % Añadimos un borde de cálulas tumorales. Como ocupan 4 casillas el
    % borde es doble. Esto va a facilitar la programación y que vamos a
    % considerar todas las casillas como si fueran interiores y no hace
    % falta distinguir los casos de las esquinas y de los bordes.
    A = [2*ones(1,n); 2*ones(1,n); A; 2*ones(1,n); 2*ones(1,n)];
    A = [2*ones(n+4,1) 2*ones(n+4,1) A 2*ones(n+4,1) 2*ones(n+4,1)];
    B = [2*ones(1,n); 2*ones(1,n); B; 2*ones(1,n); 2*ones(1,n)];
    B = [2*ones(n+4,1) 2*ones(n+4,1) B 2*ones(n+4,1) 2*ones(n+4,1)];
    m_vida = [2*ones(1,n); 2*ones(1,n); m_vida; 2*ones(1,n); 2*ones(1,n)];
    m_vida = [2*ones(n+4,1) 2*ones(n+4,1) m_vida 2*ones(n+4,1) 2*ones(n+4,1)];
    m_duplicar = [2*ones(1,n); 2*ones(1,n); m_duplicar; 2*ones(1,n); 2*ones(1,n)];
    m_duplicar = [2*ones(n+4,1) 2*ones(n+4,1) m_duplicar 2*ones(n+4,1) 2*ones(n+4,1)];
    m_contador = [2*ones(1,n); 2*ones(1,n); m_contador; 2*ones(1,n); 2*ones(1,n)];
    m_contador = [2*ones(n+4,1) 2*ones(n+4,1) m_contador 2*ones(n+4,1) 2*ones(n+4,1)];

    % Como hemos modificado el tamaño de A, hay que modificar tambien donde
    % están las células tumorales, que es en dos posiciones más de la
    % original
    for a = 1:length(fila)
        fila(a) = fila(a) + 2;
        columna(a) = columna(a) + 2;
    end
    longitud = length(fila);
    contador_vector = 0; % Se usa para añadir la posición de la nueva célula tumoral al vector de fila y columna.
    % Recorremos todas las células tumorales
    for j = 1:longitud
        % Con esto vemos si prolifera o no, la probabilidad viene dada como
        % parámetro de la función.
        prolifera = Bernu(rho);
        % Si rho=1 es que prolifera. Al igual que en mover, vemos hacia que dirección puede proliferar la célula.
        if prolifera == 1
            % Inicializamos los vectores. Estos sirven para guardar los
            % sitios donde puede proliferar la célula.
            F = zeros(1,2);
            C = zeros(1,2);
            % Contador para ver hacia cuantas posiciones se puede mover
            count=0;
            % No vamos a ver si las casillas están vacias. Lo que vamos a
            % hacer es ver donde no hay tumorales. Si pilla una CAR-T lo
            % que va a pasar es que la tumoral va a proliferar y empujará a
            % la CAR-T a otra posición.
            %1%%%%%%%%%%%%%%%%%%%%%%%%
            if A(fila(j)-2,columna(j)-2) <2 & A(fila(j)-2,columna(j)-1) <2 & A(fila(j)-1,columna(j)-2) <2 & A(fila(j)-1,columna(j)-1) <2
                    count = count + 1;
                    F(count) = fila(j)-2;
                    C(count) = columna(j)-2;
            end
            %2%%%%%%%%%%%%%%%
            if A(fila(j)-2,columna(j)-1) <2 & A(fila(j)-2,columna(j)) <2 & A(fila(j)-1,columna(j)-1) <2 & A(fila(j)-1,columna(j)) <2
                    count = count + 1;
                    F(count) = fila(j)-2;
                    C(count) = columna(j)-1;
            end
            %3%%%%%%%%%%%%%%%%%%
            if A(fila(j)-2,columna(j)) <2 & A(fila(j)-2,columna(j)+1) <2 & A(fila(j)-1,columna(j)) <2 & A(fila(j)-1,columna(j)+1) <2
                    count = count + 1;
                    F(count) = fila(j)-2;
                    C(count) = columna(j);
            end
            %4%%%%%%%%%%%%%%%%%%%%%
            if A(fila(j)-2,columna(j)+1) <2 & A(fila(j)-2,columna(j)+2) <2 & A(fila(j)-1,columna(j)+1) <2 & A(fila(j)-1,columna(j)+2) <2
                    count = count + 1;
                    F(count) = fila(j)-2;
                    C(count) = columna(j)+1;
            end
            %5%%%%%%%%%%%%%%%%%%
            if A(fila(j)-2,columna(j)+2) <2 & A(fila(j)-2,columna(j)+3) <2 & A(fila(j)-1,columna(j)+2) <2 & A(fila(j)-1,columna(j)+3) <2
                    count = count + 1;
                    F(count) = fila(j)-2;
                    C(count) = columna(j)+2;
            end
            %6%%%%%%%%%%%%%%%
             if A(fila(j)-1,columna(j)+2) <2 & A(fila(j)-1,columna(j)+3) <2 & A(fila(j),columna(j)+2) <2 & A(fila(j),columna(j)+3) <2
                    count = count + 1;
                    F(count) = fila(j)-1;
                    C(count) = columna(j)+2;
             end
             %7%%%%%%%%%%%%%%%
             if A(fila(j),columna(j)+2) <2 & A(fila(j),columna(j)+3) <2 & A(fila(j)+1,columna(j)+2) <2 & A(fila(j)+1,columna(j)+3) <2
                    count = count + 1;
                    F(count) = fila(j);
                    C(count) = columna(j)+2;
             end
             %8%%%%%%%%%%%%%%%%%%%
             if A(fila(j)+1,columna(j)+2) <2 & A(fila(j)+1,columna(j)+3) <2 & A(fila(j)+2,columna(j)+2) <2 & A(fila(j)+2,columna(j)+3) <2
                    count = count + 1;
                    F(count) = fila(j)+1;
                    C(count) = columna(j)+2;
             end
             %9%%%%%%%%%%%%%%%%%%%
             if A(fila(j)+2,columna(j)+2) <2 & A(fila(j)+2,columna(j)+3) <2 & A(fila(j)+3,columna(j)+2) <2 & A(fila(j)+3,columna(j)+3) <2
                    count = count + 1;
                    F(count) = fila(j)+2;
                    C(count) = columna(j)+2;
             end
             %10%%%%%%%%%%%%%%%%%%%
             if A(fila(j)+2,columna(j)+1) <2 & A(fila(j)+2,columna(j)+2) <2 & A(fila(j)+3,columna(j)+1) <2 & A(fila(j)+3,columna(j)+2) <2
                    count = count + 1;
                    F(count) = fila(j)+2;
                    C(count) = columna(j)+1;
             end
             %11%%%%%%%%%%%%%%%%%%%
             if A(fila(j)+2,columna(j)) <2 & A(fila(j)+2,columna(j)+1) <2 & A(fila(j)+3,columna(j)) <2 & A(fila(j)+3,columna(j)+1) <2
                    count = count + 1;
                    F(count) = fila(j)+2;
                    C(count) = columna(j);
             end
             %12%%%%%%%%%%%%%%%%%%%
             if A(fila(j)+2,columna(j)-1) <2 & A(fila(j)+2,columna(j)) <2 & A(fila(j)+3,columna(j)-1) <2 & A(fila(j)+3,columna(j)) <2
                    count = count + 1;
                    F(count) = fila(j)+2;
                    C(count) = columna(j)-1;
             end
             %13%%%%%%%%%%%%%%%%%%%
             if A(fila(j)+2,columna(j)-2) <2 & A(fila(j)+2,columna(j)-1) <2 & A(fila(j)+3,columna(j)-2) <2 & A(fila(j)+3,columna(j)-1) <2
                    count = count + 1;
                    F(count) = fila(j)+2;
                    C(count) = columna(j)-2;
             end
             %14%%%%%%%%%%%%%%%%%%%
             if A(fila(j)+1,columna(j)-2) <2 & A(fila(j)+1,columna(j)-1) <2 & A(fila(j)+2,columna(j)-2) <2 & A(fila(j)+2,columna(j)-1) <2
                    count = count + 1;
                    F(count) = fila(j)+1;
                    C(count) = columna(j)-2;
             end
             %15%%%%%%%%%%%%%%%%%%%
             if A(fila(j),columna(j)-2) <2 & A(fila(j),columna(j)-1) <2 & A(fila(j)+1,columna(j)-2) <2 & A(fila(j)+1,columna(j)-1) <2
                    count = count + 1;
                    F(count) = fila(j);
                    C(count) = columna(j)-2;
             end
             %16%%%%%%%%%%%%%%%%%%%
             if A(fila(j)-1,columna(j)-2) <2 & A(fila(j)-1,columna(j)-1) <2 & A(fila(j),columna(j)-2) <2 & A(fila(j),columna(j)-1) <2
                    count = count + 1;
                    F(count) = fila(j)-1;
                    C(count) = columna(j)-2;
             end
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             if count > 0 % Si el contador es >0 entonces la célula tiene hueco para proliferar.
                contador_vector = contador_vector + 1; % Toma el valor de 1 y se usa para añadir al vector fila y columna la posición de la nueva célula tumoral
                k = round((count-1)*rand(1,1))+1; % Generamos la posición aleatoria hacia la que se mueve
                if k > 0 
                    %1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    cont_moverT=0; % Contador de las CAR-T que hay que mover para dejar hueco a la nueva tumoral.
                    FF=zeros(1,2); % Inicializar los vectores
                    CC=zeros(1,2);
                    % Vamos a tener 4 ifs de este tipo ya que las tumorales
                    % ocupan 4 casillas, por lo que no puede haber ninguna
                    % CAR-T en ninguna de las 4 casillas que ocupa la tumoral.
                    if A(F(k),C(k)) == 1 % Esto quiere decir que en la posicion en la que va a estar la nueva tumoral hay una CAR-T que hay que mover.
                        %1%%%%%%%%%%%%%%% 
                        % Con estos ifs lo que hacemos es ver hacia que posiciones se puede mover la CAR-T
                        if A(F(k)-1,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %2%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k);
                        end
                        %3%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)+1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)+1;
                        end
                        %4%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %5%%%%%%%%%%%%%%%
                        if A(F(k),C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k);
                            CC(cont_moverT)=C(k)+2;
                        end
                        %6%%%%%%%%%%%%%%%
                        if A(F(k)+1,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+1;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %7%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %8%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)+1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)+1;
                        end
                        %9%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k);
                        end
                        %10%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %11%%%%%%%%%%%%%%%
                        if A(F(k)+1,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+1;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %12%%%%%%%%%%%%%%%
                        if A(F(k),C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k);
                            CC(cont_moverT)=C(k)-1;
                        end
                        if cont_moverT>0
                            kk=round((cont_moverT-1)*rand(1,1))+1;
                            A(FF(kk),CC(kk)) = 1;
                            B(FF(kk),CC(kk)) = B(F(k),C(k));
                            B(F(k),C(k))=0;
                            % Actualizamos las posiciones de todas las matrices que tenemos
                            m_vida(FF(kk),CC(kk)) = m_vida(F(k),C(k));
                            m_vida(F(k),C(k))=0;
                            m_duplicar(FF(kk),CC(kk)) = m_duplicar(F(k),C(k));
                            m_duplicar(F(k),C(k))=0;
                            m_contador(FF(kk),CC(kk)) = m_contador(F(k),C(k));
                            m_contador(F(k),C(k))=0;
                        end
                    end
                    %2%%%%%%%%
                    % Segunda casilla que ocupa la tumoral
                    cont_moverT=0;
                    if A(F(k),C(k)+1) == 1
                        %1%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %2%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k);
                        end
                        %3%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)+1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)+1;
                        end
                        %4%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %5%%%%%%%%%%%%%%%
                        if A(F(k),C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k);
                            CC(cont_moverT)=C(k)+2;
                        end
                        %6%%%%%%%%%%%%%%%
                        if A(F(k)+1,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+1;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %7%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %8%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)+1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)+1;
                        end
                        %9%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k);
                        end
                        %10%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %11%%%%%%%%%%%%%%%
                        if A(F(k)+1,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+1;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %12%%%%%%%%%%%%%%%
                        if A(F(k),C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k);
                            CC(cont_moverT)=C(k)-1;
                        end
                        if cont_moverT>0
                            kk=round((cont_moverT-1)*rand(1,1))+1;
                            A(FF(kk),CC(kk)) = 1;
                            B(FF(kk),CC(kk)) = B(F(k),C(k)+1);
                            B(F(k),C(k)+1)=0;
                            m_vida(FF(kk),CC(kk)) = m_vida(F(k),C(k)+1);
                            m_vida(F(k),C(k)+1)=0;
                            m_duplicar(FF(kk),CC(kk)) = m_duplicar(F(k),C(k)+1);
                            m_duplicar(F(k),C(k)+1)=0;
                            m_contador(FF(kk),CC(kk)) = m_contador(F(k),C(k)+1);
                            m_contador(F(k),C(k)+1)=0;
                        end
                    end
                    %3%%%%%%%%
                    % Tercera casilla que ocupa la tumoral
                    cont_moverT=0;
                    if A(F(k)+1,C(k)) == 1
                        %1%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %2%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k);
                        end
                        %3%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)+1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)+1;
                        end
                        %4%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %5%%%%%%%%%%%%%%%
                        if A(F(k),C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k);
                            CC(cont_moverT)=C(k)+2;
                        end
                        %6%%%%%%%%%%%%%%%
                        if A(F(k)+1,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+1;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %7%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %8%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)+1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)+1;
                        end
                        %9%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k);
                        end
                        %10%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %11%%%%%%%%%%%%%%%
                        if A(F(k)+1,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+1;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %12%%%%%%%%%%%%%%%
                        if A(F(k),C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k);
                            CC(cont_moverT)=C(k)-1;
                        end
                        if cont_moverT>0
                            kk=round((cont_moverT-1)*rand(1,1))+1;
                            A(FF(kk),CC(kk)) = 1;
                            B(FF(kk),CC(kk)) = B(F(k)+1,C(k));
                            B(F(k)+1,C(k))=0;
                            m_vida(FF(kk),CC(kk)) = m_vida(F(k)+1,C(k));
                            m_vida(F(k)+1,C(k))=0;
                            m_duplicar(FF(kk),CC(kk)) = m_duplicar(F(k)+1,C(k));
                            m_duplicar(F(k)+1,C(k))=0;
                            m_contador(FF(kk),CC(kk)) = m_contador(F(k)+1,C(k));
                            m_contador(F(k)+1,C(k))=0;
                        end
                    end
                    %4%%%%%%%%
                    % Cuarta casilla que ocupa la tumoral
                    cont_moverT=0;
                    if A(F(k)+1,C(k)+1) == 1
                        %1%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %2%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k);
                        end
                        %3%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)+1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)+1;
                        end
                        %4%%%%%%%%%%%%%%%
                        if A(F(k)-1,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)-1;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %5%%%%%%%%%%%%%%%
                        if A(F(k),C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k);
                            CC(cont_moverT)=C(k)+2;
                        end
                        %6%%%%%%%%%%%%%%%
                        if A(F(k)+1,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+1;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %7%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)+2) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)+2;
                        end
                        %8%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)+1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)+1;
                        end
                        %9%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k);
                        end
                        %10%%%%%%%%%%%%%%%
                        if A(F(k)+2,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+2;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %11%%%%%%%%%%%%%%%
                        if A(F(k)+1,C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k)+1;
                            CC(cont_moverT)=C(k)-1;
                        end
                        %12%%%%%%%%%%%%%%%
                        if A(F(k),C(k)-1) == 0
                            cont_moverT=cont_moverT+1;
                            FF(cont_moverT)=F(k);
                            CC(cont_moverT)=C(k)-1;
                        end
                        if cont_moverT>0
                            kk=round((cont_moverT-1)*rand(1,1))+1;
                            A(FF(kk),CC(kk)) = 1;
                            B(FF(kk),CC(kk)) = B(F(k)+1,C(k)+1);
                            B(F(k)+1,C(k)+1)=0;
                            m_vida(FF(kk),CC(kk)) = m_vida(F(k)+1,C(k)+1);
                            m_vida(F(k)+1,C(k)+1)=0;
                            m_duplicar(FF(kk),CC(kk)) = m_duplicar(F(k)+1,C(k)+1);
                            m_duplicar(F(k)+1,C(k)+1)=0;
                            m_contador(FF(kk),CC(kk)) = m_contador(F(k)+1,C(k)+1);
                            m_contador(F(k)+1,C(k)+1)=0;
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Añadimos la nueva célula que ha proliferado
                    A(F(k),C(k)) = 2;
                    A(F(k),C(k)+1) = 2;
                    A(F(k)+1,C(k)) = 2;
                    A(F(k)+1,C(k)+1) = 2;
                    fila(longitud+contador_vector)=F(k); % Añadimos la nueva celula tumoral al vector de tumorales
                    columna(longitud+contador_vector)=C(k);
                end
             end
        end
    end
    % Volvemos a poner las posiciones de las tumorales como al principio, es decir, sin contar la ampliación de la matriz
    for a=1:length(fila)
        fila(a)=fila(a)-2;
        columna(a)=columna(a)-2;
    end
    % Devolvemos las nuevas matrices A y B sin los bordes de seguridad que hemos añadido dentro de esta función.
    AA = A(3:end-2,3:end-2);
    BB = B(3:end-2,3:end-2);
    vida_new = m_vida(3:end-2,3:end-2);
    duplicar_new = m_duplicar(3:end-2,3:end-2);
    contador_new = m_contador(3:end-2,3:end-2);
    FILA = fila;
    COLUMNA = columna;
end