% Con esta función vamos a calcular la evolución que sigue la población de
% células CAR-T a lo largo de la simulación.
N = 100;
semilla = 1:(3000/N):3000;
% Vamos a ir actualizando una matriz en cada paso
matriz=[];
for i=1:N
    i
    [tum_final1(i),vector_CART] = fun_auto_dispersas(semilla(i),1400,20);
    matriz = [matriz; vector_CART];
end

for j=1:1200
    vector_max(j) = max(matriz((1:N),j)); %calculamos el número máximo de celulas CAR-T en ese instante de la simulación
    vector_min(j) = min(matriz((1:N),j)); %calculamos el número mínimo de celulas CAR-T en ese instante de la simulación
    vector_media(j) = mean(matriz((1:N),j)); %calculamos el número medio de células CAR-T en ese instante de la simulación
end

% Sacamos en qué momento de tiempo se encuentra el número máximo de células
% CAR-T
valor_max_cart=0;
posicion=0;
for cont=1:length(vector_max)
    if valor_max_cart < vector_max(cont)
        valor_max_cart = vector_max(cont);
        posicion=cont;
    end
end

% Creamos los vectores para pintar la región del plano correspondiente.
x1 =[0 1:1200 1200 0 1:1200 1200];
y1 =[0 vector_min 0 0 vector_max 0];

% Creamos nuevos vectores para dibujar lineas horizontales y verticales en
% el punto depico máximo.
x_horizon = 1:posicion;
y_horizon = valor_max_cart*ones(1,posicion);
x_vertical = posicion*ones(1,valor_max_cart);
y_vertical = 1:valor_max_cart;

%Representamos las tres curvas a la vez
figure(8)
% Relleno entre los valores máx y mín
fill(x1,y1,'cyan','FaceAlpha',0.3)
hold on
% Curva de valores máximos
plot(1:1200,vector_max, 'r','LineWidth',2)
% Curva de valores mínimos
plot(1:1200,vector_min, 'g','LineWidth',2)
% Curva de valores medeios
plot(1:1200,vector_media, 'b:','LineWidth',2)
% Punto para representar el pico máximo
plot(posicion,valor_max_cart, 'k*','LineWidth',2)
% Líneas horizontales y verticales
plot(x_horizon,y_horizon,'k-.')
plot(x_vertical,y_vertical,'k-.')
title('Evolución del número de células CAR-T')
subtitle('CAR-T_0 = 20, Ntum_0 = 1400')
legend('','Nº máx CAR-T', 'Nº mín CAR-T', 'Nº medio CAR-T','Pico máximo CAR-T')
xlabel('Tiempo (horas)') 
ylabel('Número de células CAR-T') 
hold off


