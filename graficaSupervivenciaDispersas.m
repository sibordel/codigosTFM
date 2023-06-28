% vectores que guarda los datos de los tumores que crecen
vX1 = [-1];
vY1 = [-1];
% vectores que guarda los datos de los tumores que se reducen entre un 0 y un 25 % 
vX2 = [-1];
vY2 = [-1];
% vectores que guarda los datos de los tumores que se reducen entre un 25 y un 50 % 
vX3 = [-1];
vY3 = [-1];
% vectores que guarda los datos de los tumores que se reducen entre un 50 y un 75 % 
vX4 = [-1];
vY4 = [-1];
% vectores que guarda los datos de los tumores que se reducen entre un 75 y un 100 % 
vX5 = [-1];
vY5 = [-1];
% vectores que guarda los datos de los tumores que desaparecen % 
vX6 = [-1];
vY6 = [-1];

N = 20;
semilla = 1:(3000/N):3000;
% Para que ejecute más rápido en el ordenador y haga las iteraciones en paralelo
%parpool(16)
parfor contCART=1:100
    for contTUM=0:32
        ntum0 = 600 + contTUM*25;
        vector_tum_final = [];
        for iter=1:N % Calculamos 20 iteraciones por cada caso
            [vector_tum_final(iter),vector_CART] = fun_auto_dispersas(semilla(iter),ntum0,contCART);
            fprintf('contCART: %d, contTUM: %d, iter: %d\n',contCART,contTUM,iter)
        end
        % Guardamos en un vector el número final de tumorales en cada una de las 100 iteraciones
        tum_final=mean(vector_tum_final); % Calculamos la media de tumorales finales para ver los porcentajes en los que se ha reducido o aumentado
        if tum_final > ntum0 % En este caso el tumor aumenta
            vX1 = [vX1 contCART];
            vY1 = [vY1 ntum0];
        elseif tum_final <= ntum0 & tum_final > 0.75*ntum0 % En este caso el tumor se reduce entre un 0 y un 25%
            vX2 = [vX2 contCART];
            vY2 = [vY2 ntum0];
            elseif tum_final <= 0.75*ntum0 & tum_final > 0.5*ntum0 % En este caso el tumor se reduce entre un 25 y un 50%
            vX3 = [vX3 contCART];
            vY3 = [vY3 ntum0];
            elseif tum_final <= 0.5*ntum0 & tum_final > 0.25*ntum0 % En este caso el tumor se reduce entre un 50 y un 75%
            vX4 = [vX4 contCART];
            vY4 = [vY4 ntum0];
            elseif tum_final <= 0.75*ntum0 & tum_final > 0 % En este caso el tumor se reduce entre un 75 y un 99%
            vX5 = [vX5 contCART];
            vY5 = [vY5 ntum0];
            elseif tum_final == 0 % En este caso el tumor desaparece
            vX6 = [vX6 contCART];
            vY6 = [vY6 ntum0];
        end
    end
end

% % save vX1;
% % save vY1;
% % save vX2;
% % save vY2;
% % save vX3;
% % save vY3;
% % save vX4;
% % save vY4;
% % save vX5;
% % save vY5;
% % save vX6;
% % save vY6;

%Dibujamos la gráfica
figure (1)
if length(vX1) > 0 
scatter(vX1,vY1,'r','filled')
end

hold on

if length(vX2) > 0 
    scatter(vX2,vY2,'m','filled')
end
if length(vX3) > 0 
    scatter(vX3,vY3,'y','filled')
end
if length(vX4) > 0 
    scatter(vX4,vY4,'c','filled')
end
if length(vX5) > 0 
    scatter(vX5,vY5,'b','filled')
end
if length(vX6) > 0 
    scatter(vX6,vY6,'g','filled')
end
title('Supervivencia del tumor en función del número inicial de células (CAR-T y tumorales)')
xlabel('Nº inicial de células CAR-T') 
ylabel('Nº inicial de células tumorales')
xlim([0 105])
ylim([500 1500])
legend('El tumor crece','El tumor se reduce menos de un 25%', 'El tumor se reduce entre un 25-50%','El tumor se reduce entre un 50-75%','El tumor se reduce más de un 75%','El tumor desaparece')
hold off
