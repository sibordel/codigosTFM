tic
n = 200; % nxn tamaño de la matriz
ntum = 1000; % nº de celulas infectadas
A = zeros(n+1,n+1); % guardamos los tipos de células: tumoral o CAR-T
B = zeros(n+1,n+1); % guardamos la carga de las CAR-T
m_vida = zeros(n+1,n+1); % guardamos el tiempo de vida de cada CAR-T
m_duplicar = zeros(n+1,n+1); % guardamos qué células están en fase de duplicación
m_contador = zeros(n+1,n+1); % guardamos el número de duplicaciones que lleva cada CAR-T
k=0; % contador de celulas tumorales

% Para que la simulación siempre sea la misma fijamos una semilla de aleatoriedad
rng('default');
rng(2000); % el número es la semilla


fila=ceil(n/2);
col=ceil(n/2);


% En fila1, col1 se guardan las posiciones donde se encuentran las células
% tumorales. Como cada tumoral ocupa 4 casillas vamos a guardar la primera
% posición, es decir, la casilla de arriba a la izquierda. A partir de esta
% posición podemos sacar las demás sin problema.

% Las matrices son de tamaño (n+1)x(n+1) pero las células solo están entre
% las posiciones 1 y n. Esto hace que se forme un borde en el lado derecho
% y abajo, es decir, no puede haber cells en la última fila ni en la última
% columna. Luego añadiremos el borde superior y el borde inferior.

col1 = 1:ntum;
fila1 = 1:ntum;
while k < ntum
    fila = ceil((n-1)*rand(1,1)); % Generamos un posición de la matriz
    col = ceil((n-1)*rand(1,1));
    
    % Rellenar con células infectadas
    if A(fila,col) == 0 % si el valor de la posición es 0, añadimos la célula infectada
        A(fila, col) = 1;
        A(fila+1, col) = 1;
        A(fila, col+1) = 1;
        A(fila+1, col+1) = 1;
        fila1(k+1) = fila;
        col1(k+1) = col;
        if fila>1 & col>1   % Rellenamos lo de alrededor de -1s para saber que no puede caer ahí la próxima célula
            if A(fila-1,col-1) == 0
                A(fila-1,col-1) = -1;
            end
        end
        if fila>1
            if A(fila-1, col) == 0
               A(fila-1, col) = -1;
            end
            if A(fila-1, col+1) == 0
               A(fila-1, col+1) = -1;
            end
        end
        if col>1
            if A(fila, col-1) == 0
               A(fila, col-1) = -1;
            end
            if A(fila+1, col-1) == 0
               A(fila+1, col-1) = -1;
            end
        end
    end
    k = (numel(A(A==1)))/4; % Contar nº de células que llevamos insertadas
end

A = A + abs(A); % Quitamos los -1s (Las células infectadas ahora son 2s)


% Rellenar con células CAR-T
nT = 20;
kT=0; % contador de celulas CAR-T
while kT < nT
    filaT = ceil(n*rand(1,1)); % Generamos un posición de la matriz
    colT = ceil(n*rand(1,1));
    if A(filaT,colT) == 0 % si el valor de la posición es 0, añadimos la celula CAR-T
        A(filaT, colT) = 1;
        B(filaT,colT) = 9;
        m_vida(filaT,colT) = 336;
    end  
    kT = numel(A(A==1)); % Contar nº de células que llevamos insertadas
end

% Borde de seguridad para no salirnos de rango en las matrices
% Esto me añade una fila de ceros al principio
A = [zeros(1,n+1); A];
% Esto me añade una columna de ceros al principio
A = [zeros(n+2,1) A];
% Esto me añade una fila de ceros al principio
B = [zeros(1,n+1); B];
% Esto me añade una columna de ceros al principio
B = [zeros(n+2,1) B];
% Añadimos el borde de seguridad en el resto de matrices
m_vida = [zeros(1,n+1); m_vida];
m_vida = [zeros(n+2,1) m_vida];
m_duplicar = [zeros(1,n+1); m_duplicar];
m_duplicar = [zeros(n+2,1) m_duplicar];
m_contador = [zeros(1,n+1); m_contador];
m_contador = [zeros(n+2,1) m_contador];

% Graficamos los movimientos de las células CAR-T, disparar, matar, mover y proliferar células tumorales
T =1200; % tiempo en horas
for tiempo = 1:T
      % Carga citotoxica acumulada en cada célula tumoral. Vector fila que
      % tiene tantas componentes como células tumorales el autómata.
      % Todas sus componentes empiezan inicializadas en 0.
      carga_toxica_ctumoral = zeros(1,length(fila1));
      % Para definir el color del autómata.
      COLOR=[178/255 218/255 250/255; 147/255 112/255 219/255; 0.6350 0.0780 0.1840];
      clim([0 2])
      colormap(COLOR)
      imagesc(A(2:end-1,2:end-1))
      title(['Hora ' num2str(tiempo) ' de ' num2str(T)])

      % Actualizamos el valor de carga citotóxica que tiene cada célula tumoral.
      % Hemos añadido la m_duplicar a disparar. Aqui sabemos que células
      % CAR-T están en fase de duplicación y no pueden moverse
      [carga_toxica_ctumoral,B,m_duplicar,m_contador] = disparar(fila1, col1, carga_toxica_ctumoral, B, m_duplicar,m_contador, 24,20);
      fila_new = [];
      col_new = [];
      carga_toxica_ctumoral_new = [];

      % Vemos si la célula tumoral muere o no muere
      for i = 1:length(carga_toxica_ctumoral)
        % La función morir dependiendo de la carga citotoxica que haya
        % recibido cada célula tumoral asigna distintas probabilidades de
        % muerte. Cuanta más carga tenga, más probabilidad de morir. Le
        % pasamos a morir célula tumoral por célula tumoral
        matamos = morir(carga_toxica_ctumoral(i));
        if matamos == 1
            % Este es el caso en el que la tumoral muere, por lo que
            % para eliminarla ponemos su posición a 0. Hay que tener
            % cuidado ya que la matriz A está ampliada, por eso tenemos que
            % sumar un +1 más al correspondiente.
            A(fila1(i)+1, col1(i)+1) = 0;
            A(fila1(i)+2, col1(i)+1) = 0;
            A(fila1(i)+1, col1(i)+2) = 0;
            A(fila1(i)+2, col1(i)+2) = 0;
        else
            % En este caso la célula no muere por lo que actualizamos el
            % valor de fila_new, col_new y carga_new. Solo nos guardamos la
            % posición y la carga de las tumorales restantes.
            fila_new = [fila_new fila1(i)];
            col_new = [col_new col1(i)];
            carga_toxica_ctumoral_new = [carga_toxica_ctumoral_new carga_toxica_ctumoral(i)];
        end
      end
      % Una vez hemos actualizado todas las muertes de las tumorales lo que
      % hacemos es actualizar también los vectores fila1, col1 y
      % carga_toxica_ctumoral para usarlos en la siguiente iteración del bucle.
      fila1 = fila_new;
      col1 = col_new;
      carga_toxica_ctumoral = carga_toxica_ctumoral_new;

        
      % Una vez hayan disparado y matado, las células CAR-T se mueven por el
      % autómata. Cuidado con las dimensiones de la matriz ya que tenemos A
      % y B ampliadas y a mover le tenemos que pasar matrices nxn.
      [A,B,m_vida,m_duplicar,m_contador] = mover(A(2:end-1,2:end-1),B(2:end-1,2:end-1),m_vida(2:end-1,2:end-1),m_duplicar(2:end-1,2:end-1),m_contador(2:end-1,2:end-1));
      
      % Recargamos la CAR-T
      [B] = recargar(A,B,m_duplicar);

      % Una vez se hayan movido las células CAR-T, las tumorales proliferan
      [A,B,m_vida,m_duplicar,m_contador,fila1,col1] = proliferar(fila1,col1,A,B,m_vida,m_duplicar,m_contador,0.001);
      % Volvemos a definir el borde de seguridad.
      A = [zeros(1,n); A; zeros(1,n)];
      A = [zeros(n+2,1) A zeros(n+2,1)];
      B = [zeros(1,n); B; zeros(1,n)];
      B = [zeros(n+2,1) B zeros(n+2,1)];
      m_vida = [zeros(1,n); m_vida; zeros(1,n)];
      m_vida = [zeros(n+2,1) m_vida zeros(n+2,1)];
      m_duplicar = [zeros(1,n); m_duplicar; zeros(1,n)];
      m_duplicar = [zeros(n+2,1) m_duplicar zeros(n+2,1)];
      m_contador = [zeros(1,n); m_contador; zeros(1,n)];
      m_contador = [zeros(n+2,1) m_contador zeros(n+2,1)];
      
      % Vamos a ver cuando las CAR-T de duplican
      % Primero vemos que valores son igual a 1. Estas cells son las que se van a duplicar
      [A,B,m_vida,m_duplicar,m_contador] = duplicar(A(2:end-1,2:end-1),B(2:end-1,2:end-1),m_vida(2:end-1,2:end-1),m_duplicar(2:end-1,2:end-1),m_contador(2:end-1,2:end-1));
      % Volvemos a definir el borde de seguridad.
      A = [zeros(1,n); A; zeros(1,n)];
      A = [zeros(n+2,1) A zeros(n+2,1)];
      B = [zeros(1,n); B; zeros(1,n)];
      B = [zeros(n+2,1) B zeros(n+2,1)];
      m_vida = [zeros(1,n); m_vida; zeros(1,n)];
      m_vida = [zeros(n+2,1) m_vida zeros(n+2,1)];
      m_duplicar = [zeros(1,n); m_duplicar; zeros(1,n)];
      m_duplicar = [zeros(n+2,1) m_duplicar zeros(n+2,1)];
      m_contador = [zeros(1,n); m_contador; zeros(1,n)];
      m_contador = [zeros(n+2,1) m_contador zeros(n+2,1)];

      % Restamos una hora a la fase de duplicación de las demás CAR-T
      [f_dup2,c_dup2] = find(m_duplicar > 1);
      lon2 = length(f_dup2);
      if lon2 > 0
          for uu=1:lon2
              m_duplicar(f_dup2(uu),c_dup2(uu)) = m_duplicar(f_dup2(uu),c_dup2(uu)) -1;
          end
      end

      % Primero tenemos que eliminar los que sean igual a 1 y luego
      % restarles una hora de vida a los demás.
      % Los que sean igual a 1 los elimamos de la matriz. Esto quiere decir
      % que en el siguiente paso ya no hay CAR-T
      [f_muerte,c_muerte] = find(m_vida==1);
      longitud = length(f_muerte);
      if longitud > 0
          for r=1:longitud
              m_vida(f_muerte(r),c_muerte(r))=0;
              A(f_muerte(r),c_muerte(r))=0;
              B(f_muerte(r),c_muerte(r))=0;
              m_duplicar(f_muerte(r),c_muerte(r))=0;
          end
      end

      % Vamos a restarles una hora de vida
      [f_vida,c_vida]=find(m_vida > 1);
      longi = length(f_vida);
      if longi > 0
          for p = 1:longi
              m_vida(f_vida(p),c_vida(p)) = m_vida(f_vida(p),c_vida(p)) - 1;
          end
      end
      
      drawnow
      % Para que la simulación vaya más lenta
      %pause(0.2)
      % Contamos la cantidad de tumorales en cada tiempo para hacer un
      % gráfico de cómo van disminuyendo a medida que pasa la simulación.
      evolucion(tiempo) = length(carga_toxica_ctumoral);
      num_cell_t(tiempo) = numel(A(A==1));
end

      % Sacamos en qué momento de tiempo se encuentra el número máximo de células tumorales
      valor_max_tum=0;
      posicionTUM=0;
      for cont=1:length(evolucion)
          if valor_max_tum < evolucion(cont)
            valor_max_tum = evolucion(cont);
            posicionTUM=cont;
          end
      end

      % Creamos nuevos vectores para dibujar líneas horizontales y verticales en
      % el punto de pico máximo.
      x_horizonTUM = 1:posicionTUM;
      y_horizonTUM = valor_max_tum*ones(1,posicionTUM);
      x_verticalTUM = posicionTUM*ones(1,valor_max_tum);
      y_verticalTUM = 1:valor_max_tum;

      % Gráfica de cómo van disminuyendo las tumorales.
      figure(2)
      plot(0:T, [ntum evolucion],'m*-', 'LineWidth',2)
      hold on
      % Punto para representar el pico máximo
      plot(posicionTUM,valor_max_tum, 'k*','LineWidth',2)
      % Líneas horizontales y verticales
      plot(x_horizonTUM,y_horizonTUM,'k-.')
      plot(x_verticalTUM,y_verticalTUM,'k-.')
      title("Evolución células tumorales");
      subtitle('CAR-T_0 = 20, Ntum_0 = 1000')
      xlabel('Tiempo (horas)') 
      ylabel('Nº de células tumorales')
      xlim([0 T])
      hold off


      % Sacamos en qué momento de tiempo se encuentra el número máximo de CAR-T
      valor_max_cart=0;
      posicionCART=0;
      for cont=1:length(num_cell_t)
          if valor_max_cart < num_cell_t(cont)
            valor_max_cart = num_cell_t(cont);
            posicionCART=cont;
          end
      end

      % Creamos nuevos vectores para dibujar líneas horizontales y verticales en
      % el punto de pico máximo.
      x_horizonCART = 1:posicionCART;
      y_horizonCART = valor_max_cart*ones(1,posicionCART);
      x_verticalCART = posicionCART*ones(1,valor_max_cart);
      y_verticalCART = 1:valor_max_cart;


      %Grafica de las CAR-T
      figure(3)
      plot(0:T, [nT num_cell_t],'b*-', 'LineWidth',2)
      hold on
      % Punto para representar el pico máximo
      plot(posicionCART,valor_max_cart, 'k*','LineWidth',2)
      % Líneas horizontales y verticales
      plot(x_horizonCART,y_horizonCART,'k-.')
      plot(x_verticalCART,y_verticalCART,'k-.')
      title("Evolución CAR-T");
      subtitle('CAR-T_0 = 20, Ntum_0 = 1000')
      xlabel('Tiempo (horas)') 
      ylabel('Nº de células tumorales')
      xlim([0 T])
      hold off
      toc

      % Esto es para usar la función bwconncomp y contar el número de
      % colonias que quedan de células tumorales
      % No tenemos una matriz que tenga solo las células tumorales así que
      % la vamos a crear. Vamos a crear una nueva matriz TUM que solo tenga el número
      % de células tumorales finales. Para ello, vamos a restar 1 a las
      % posiciones de la matriz A que sean CAR-T cells (es decir las
      % posiciones en las que hay 1). Las posiciones de las tumorales son 2
%       TUM = A;
%       [f_cart,c_cart] = find(TUM == 1);
%       longi3 = length(f_cart);
%       if longi3 > 0
%           for pp = 1:longi3
%               TUM(f_cart(pp),c_cart(pp)) = 0;
%           end
%       end
%       bwconncomp(TUM)
