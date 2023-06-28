N = 100;
semilla = 1:(3000/N):3000;
%parpool(16)
parfor i=1:N
    i
    [tum_final1(i),NumCarT] = fun_auto_dispersas(semilla(i),600,20);
    [tum_final2(i),NumCarT] = fun_auto_dispersas(semilla(i),800,20);
    [tum_final3(i),NumCarT] = fun_auto_dispersas(semilla(i),1000,20);
    [tum_final4(i),NumCarT] = fun_auto_dispersas(semilla(i),1200,20);
    [tum_final5(i),NumCarT] = fun_auto_dispersas(semilla(i),1400,20);
end
% Calculamos las medias de todos los casos.
media_normal1 = mean(tum_final1);
media_normal2 = mean(tum_final2);
media_normal3 = mean(tum_final3);
media_normal4 = mean(tum_final4);
media_normal5 = mean(tum_final5);
% Representación gráfica
figure(6)
scatter(1:N,tum_final1, 'c','filled')
hold on
scatter(1:N,tum_final2, 'm','filled')
scatter(1:N,tum_final3, 'g','filled')
scatter(1:N,tum_final4, 'r','filled')
scatter(1:N,tum_final5, 'y','filled')
% Representamos líneas horizontales con las medias
yline(media_normal1,':c',LineWidth=2)
yline(media_normal2,':m',LineWidth=2)
yline(media_normal3,':g',LineWidth=2)
yline(media_normal4,':r',LineWidth=2)
yline(media_normal5,':y',LineWidth=2)
title('Número de células tumorales finales')
subtitle('CAR-T_0 = 20')
legend('ntum_0 = 600','ntum_0 = 800','ntum_0 = 1000','ntum_0 = 1200','ntum_0 = 1400')
xlabel('Nº simulaciones') 
ylabel('Nº final de células tumorales')
hold off
