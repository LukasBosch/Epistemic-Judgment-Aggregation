function w = calcScenarioWeight(priorProb, Facts)

%   Calculates scenarioweights for given prior probabilities of premises
%			
%   Code from the Bachelor Thesis: "Quantitative Epistemologische Demokratie:
% Simulation und Vergleich der konklusions- und der Prämissenbasierten
% Judgment Aggregation", Lukas Bosch, KIT



    prems = length(Facts)-1; %prems = premises
    prob = priorProb;
    for j = 1:prems
        if Facts(j) == 0
            prob(j) = 1 - priorProb(j);
        end
    end
    w = 1;
    for j = 1:prems
        w = w * prob(j);
    end
    w;
end  

