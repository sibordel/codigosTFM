function [pium,BB,m_duplicar_new,contador_new] = disparar(fila1, columna1, cargatoxicactumoral, B, m_duplicar,m_contador,tiempo_dup,dup_max)

% El parámetro tiempo_prof hace referencia al tiempo que tardan las CAR-T
% en duplicarse una vez que hayan disparado.

% Hemos metido B ampliada entonces las posiciones de las tumorales cambian y
% las de las células CAR-T también. Porque la fila1 y columan1 que metemos están
% calculadas con el tamaño inicial de la matriz.

% Meter m_duplicar en las mismas condiciones que B

% Con este bucle vamos a recorrer todas las tumorales. Como una tumoral
% ocupa 4 casillas, puede estar en contacto con 12 posibles células CAR-T.
% Vamos a ir recorriendo estas 12 posiciones posibles y sumarle la carga
% que tengan. Si en alguna posición no hay CAR-T entonces no se sumará
% nada. Si hay CAR-T se sumará la carga que tenga disponible.

% La carga de las CAR-T tiene 8 niveles, para que no hubiera problemas de
% programación esta carga está definida en el intervalo del 1 al 9. Por lo
% tanto, cuando le sumemos la carga a la tumoral tenemos que restarle 1
% para trasladarlo al intervalo [0,8].

% Si entra en el if es porque hay celula CAR-T y se ha producido un disparo. Si
% no entra en el if es porque o no hay nada o la CAR-T está vaciada.

% Vamos a considerar que las células CAR-T disparan cuando tienen la carga
% al máximo. Si queremos cambiar este comportamiento y que disparen siempre
% que puedan aunque no sea con la carga máxima, solo tendríamos que cambiar
% el valor de la variable carga_con_la_que_disparan.
carga_con_la_que_disparan = 8;

    for i = 1:length(fila1)
        if B(fila1(i),columna1(i)) > carga_con_la_que_disparan
            cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila1(i),columna1(i))-1;
            if m_contador(fila1(i),columna1(i)) < dup_max
                m_duplicar(fila1(i),columna1(i)) = tiempo_dup;
            end
        end
        if B(fila1(i),columna1(i)+1) > carga_con_la_que_disparan
            cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila1(i),columna1(i)+1)-1;
            if m_contador(fila1(i),columna1(i)+1) < dup_max
                m_duplicar(fila1(i),columna1(i)+1) = tiempo_dup;
            end
        end
        if B(fila1(i),columna1(i)+2) > carga_con_la_que_disparan
            cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila1(i),columna1(i)+2)-1;
            if m_contador(fila1(i),columna1(i)+2) < dup_max
                m_duplicar(fila1(i),columna1(i)+2) = tiempo_dup;
            end
        end
        if B(fila1(i),columna1(i)+3) > carga_con_la_que_disparan
            cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila1(i),columna1(i)+3)-1;
            if m_contador(fila1(i),columna1(i)+3) < dup_max
                m_duplicar(fila1(i),columna1(i)+3) = tiempo_dup;
            end
        end
        %%%%%%%%%%%%%%%%
        if B(fila1(i)+1,columna1(i)) > carga_con_la_que_disparan
            cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila1(i)+1,columna1(i))-1;
            if m_contador(fila1(i)+1,columna1(i)) < dup_max
                m_duplicar(fila1(i)+1,columna1(i)) = tiempo_dup;
            end
        end
        if B(fila1(i)+1,columna1(i)+3) > carga_con_la_que_disparan
            cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila1(i)+1,columna1(i)+3)-1;
            if m_contador(fila1(i)+1,columna1(i)+3) < dup_max
                m_duplicar(fila1(i)+1,columna1(i)+3) = tiempo_dup;
            end
        end
        %%%%%%%%%%%%%%%%
        if B(fila1(i)+2,columna1(i)) > carga_con_la_que_disparan
            cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila1(i)+2,columna1(i))-1;
            if m_contador(fila1(i)+2,columna1(i)) < dup_max
                m_duplicar(fila1(i)+2,columna1(i))= tiempo_dup;
            end
        end
        if B(fila1(i)+2,columna1(i)+3) > carga_con_la_que_disparan
            cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila1(i)+2,columna1(i)+3)-1;
            if m_contador(fila1(i)+2,columna1(i)+3) < dup_max
                m_duplicar(fila1(i)+2,columna1(i)+3)= tiempo_dup;
            end
        end
        %%%%%%%%%%%%%%%
        if B(fila1(i)+3,columna1(i)) > carga_con_la_que_disparan
            cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila1(i)+3,columna1(i))-1;
            if m_contador(fila1(i)+3,columna1(i)) < dup_max
                m_duplicar(fila1(i)+3,columna1(i))= tiempo_dup;
            end
        end
        if B(fila1(i)+3,columna1(i)+1) > carga_con_la_que_disparan
            cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila1(i)+3,columna1(i)+1)-1;
            if m_contador(fila1(i)+3,columna1(i)+1) < dup_max
                m_duplicar(fila1(i)+3,columna1(i)+1)= tiempo_dup;
            end
        end
        if B(fila1(i)+3,columna1(i)+2) > carga_con_la_que_disparan
            cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila1(i)+3,columna1(i)+2)-1;
            if m_contador(fila1(i)+3,columna1(i)+2) < dup_max
                m_duplicar(fila1(i)+3,columna1(i)+2)= tiempo_dup;
            end
        end
        if B(fila1(i)+3,columna1(i)+3) > carga_con_la_que_disparan
            cargatoxicactumoral(i)= cargatoxicactumoral(i)+B(fila1(i)+3,columna1(i)+3)-1;
            if m_contador(fila1(i)+3,columna1(i)+3) < dup_max
                m_duplicar(fila1(i)+3,columna1(i)+3)= tiempo_dup;
            end
        end

        % Ponemos a 0 la carga de las células CAR-T. Como han disparado, se quedan sin carga.
        B(fila1(i),columna1(i))=0;
        B(fila1(i),columna1(i)+1)=0;
        B(fila1(i),columna1(i)+2)=0;
        B(fila1(i),columna1(i)+3)=0;
        B(fila1(i)+1,columna1(i))=0;
        B(fila1(i)+1,columna1(i)+3)=0;
        B(fila1(i)+2,columna1(i))=0;
        B(fila1(i)+2,columna1(i)+3)=0;
        B(fila1(i)+3,columna1(i))=0;
        B(fila1(i)+3,columna1(i)+1)=0;
        B(fila1(i)+3,columna1(i)+2)=0;
        B(fila1(i)+3,columna1(i)+3)=0;
    end
    % Devolvemos los valores calculados
    pium = cargatoxicactumoral;
    BB = B;
    m_duplicar_new = m_duplicar;
    contador_new = m_contador;
end