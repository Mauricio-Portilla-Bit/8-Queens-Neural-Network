clear, clc


%SOLUCIÓN EN MATLAB DEL PROBLEMA DE LAS 8 REINAS  

%Mauricio Portilla Ramírez     | A01284377   


% TIEMPO ESTIMADO: 0.003 segundos 


tic 

% 1) Inicialización
Niter = 10000; 
NPopulation = 100; 
NBoard = 8;
Pmut = 1; % Probabilidad de Mutación 
Score = zeros(NPopulation,1);
Punto_de_cruzamiento = 0; 
Padres = zeros(2, NBoard);
minJaques = [100,100]; 
maxJaques = [0,0]; 
numHijos = 10; 
hijos = zeros(numHijos, NBoard);
mejoresHijos = zeros(2, NBoard); 
posMutacion = 0; 
ScoreHijos = zeros(numHijos, 1);
minJaquesHijos = [100,100]; 
posDebiles = [0,0];

% 2) Generación de Población Inicial 
Population = zeros(NPopulation,NBoard); 

for i=1:NPopulation
    Population(i,1:NBoard) = randperm(NBoard);
end


% 3) Evaluar función de aptitud 
for k=1:NPopulation
    Score(k) = contarJaques(Population(k,1:NBoard), NBoard);
    if Score(k) == 0 
        break; 
    end 
end

% 4) Ciclo de Algorítmo Genético

if Score(k) ~= 0
    
    for p=1:Niter
        
        % Generar Población 
        for i=1:NPopulation
            Population(i,1:NBoard) = randperm(NBoard);  
        end

        % Selección de los más aptos 
        for k=1:NPopulation
            Score(k) = contarJaques(Population(k,1:NBoard), NBoard);
            if Score(k) == 0 
                break; 
            end
            % Los Cromsomas más fuertes 
            if Score(k) < minJaques(1) || Score(k) < minJaques(2)
                if minJaques(1) > minJaques(2)
                    Padres(1,1:NBoard) = Population(k,1:NBoard); 
                    minJaques(1) = Score(k);
                else
                    Padres(2,1:NBoard) = Population(k,1:NBoard); 
                    minJaques(2) = Score(k);
                end
            end
            
            % Los Cromosomas más débiles 
            if Score(k) > maxJaques(1) || Score(k) > maxJaques(2)
                if maxJaques(1) < maxJaques(2)
                    maxJaques(1) = Score(k);
                    posDebiles(1) = k;
                else
                    maxJaques(2) = Score(k);
                    posDebiles(2) = k; 
                end
            end            
            
        end
        
        if Score(k) == 0 
            break; 
        end
        
        % Reproducción 
        for r=1:numHijos
            Punto_de_cruzamiento = randi(NBoard);
            hijos(r,1:Punto_de_cruzamiento) = Padres(1,1:Punto_de_cruzamiento);
            hijos(r,Punto_de_cruzamiento:NBoard) = Padres(2,Punto_de_cruzamiento:NBoard);
        end
        
        % Mutación 
        for r=1:numHijos
            if (randi(10) > 2)
                posMutacion = randi(NBoard-1); 
                dummyP = hijos(r,posMutacion);
                hijos(r, posMutacion) = hijos(r, posMutacion + 1); 
                hijos(r, posMutacion + 1) = dummyP;
            end
        end
           
        % Evaluación de Función 
        for k=1:numHijos
            ScoreHijos(k) = contarJaques(hijos(k,1:NBoard), NBoard);
            if Score(k) == 0 
                break; 
            end
            
            if ScoreHijos(k) < minJaquesHijos(2) || ScoreHijos(k) < minJaquesHijos(1)
                if minJaquesHijos(2) < minJaquesHijos(1)
                    mejoresHijos(1,1:NBoard) = hijos(k,1:NBoard);
                else
                    mejoresHijos(2,1:NBoard) = hijos(k,1:NBoard);
                end
            end
            
        end
           
        if Score(k) == 0 
            break; 
        end
        
        % Eliminar a los individuos menos aptos 
         Population(posDebiles(1), 1:NBoard) = mejoresHijos(1,1:NBoard);
         Population(posDebiles(2), 1:NBoard) = mejoresHijos(2,1:NBoard);
    end

end
toc

% Plot
if Score(k) == 0 
    Population(k,1:NBoard)
    scatter([1,2,3,4,5,6,7,8],Population(k,1:NBoard))
    title("Tablero de Ajedrez")
    xlabel("X")
    ylabel("Y")
    grid on 
end


    
% 3) Evaluar la función de Aptitud 
function [Jaques]  = contarJaques(V, NBoard)

    Jaques =  0;
    % Por la naturaleza del genoma, no existen elementos en la misma columa
    % (e inicialmente, tampoco en la misma fila), por lo que no existen
    % reinas en el mismo punto bajo ninguna circunstancia. 
    
     for i=1:NBoard
        
        % Evaluar si hay una reina en la misma fila
        for j=i:NBoard-i
            if V(i)== V(j + 1)
                Jaques = Jaques + 1;  
            end
        end
        
        % Evaluar si hay un elemento en la misma diagonal (de abajo a arriba)
        for j=1:NBoard - V(i)
            
            % Si se sale por la derecha
            if i + j > NBoard
               break;  
            end
            
            % Si existe un elemento en la diagonal, contar jaque
            if V(i) + j == V(i + j)
                Jaques = Jaques + 1;  
            end
        end
        
        % Evaluar si hay un elemento en la misma diagonal (de arriba a abajo)
        for j=1:V(i) - 1 
            
            % Si se sale por la derecha  
            if i + j > NBoard
               break; 
            end
            
            % Si existe un elemento en la diagonal, contar jaque
            if V(i) - j == V(i + j)
                Jaques = Jaques + 1;  
            end
            
        end
        
     end
        
    % Evaluar si hay un elemento en la misma diagonal (de derecha a izquierda)
end




