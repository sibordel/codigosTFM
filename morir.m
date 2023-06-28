function muerto = morir(cargatoxicactumoral)
% Le llega la carga citotóxica que tiene cada célula tumoral y aplica una
% Bernoulli en función de las probabilidades de muerte asignadas.
    muerto = 0;
    if cargatoxicactumoral> 0 & cargatoxicactumoral <9
        muerto = Bernu(0.05);
    end
    if cargatoxicactumoral> 8 & cargatoxicactumoral <17
        muerto = Bernu(0.12);
    end
    if cargatoxicactumoral> 16 & cargatoxicactumoral < 33
        muerto = Bernu(0.5);
    end
    if cargatoxicactumoral> 32 & cargatoxicactumoral < 49
        muerto = Bernu(0.8);
    end
    if cargatoxicactumoral> 48
        muerto = Bernu(0.99);
    end    
end


