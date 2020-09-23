function m = majority(V, voters) 

%   Calculates the majority votes (Judgment-Sets on the premises)
%			
%   Code from the Bachelor Thesis: "Quantitative Epistemologische Demokratie:
% Simulation und Vergleich der konklusions- und der Prämissenbasierten
% Judgment Aggregation", Lukas Bosch, KIT




    %V steht für Votetable. gemeint ist Matrix mit allen Voter-votes die
    %ausgewertet werden sollen
   	m = sum(V)/voters;
    m = m > 0.5;
    if voters == 1
        m = V;
    end
end

