function [SMatrix, string_SMatrix] = genScens(premises,resultRelation)

%   Calculates truthtable for premises and given conclusion and the
%   truthtable in String-form for GUI scneario dropdown
%			
%   Code from the Bachelor Thesis: "Quantitative Epistemologische Demokratie:
% Simulation und Vergleich der konklusions- und der Prämissenbasierten
% Judgment Aggregation", Lukas Bosch, KIT

    scens = 2^premises;
    %Deklariert Symbolische Variablen    
    syms A B C D E F G H I J K L M N O;
    %A&B, definiert symbolischen Ausdruck aus String rR
    conclusion(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O) = str2sym(resultRelation);
    %Übersetzt symbolischen Ausdruck in nutzbare Matlabfunktion
	conclude = matlabFunction(conclusion, "Vars",{[A B C D E F G H I J K L M N O]});
            
    %Jetzt gehts erst "wirklich" los
    %Generiert Wahrheitstafel (Alle möglichen Belegungen der Prämissen)
    SMatrix = (dec2bin(0:((2^premises)-1),premises)=='1');     
    %S = logical(S);
    %ScenarioConclusions, Vektor, dessen i-tes Element die Conclusion der i-ten Zeile von S ist.
    SC = conclude(SMatrix)
    %Konvertiert nach logical und hängt Vektor SCo an Matrix S an.
    SC = logical(SC);           
    %Hängt Spaltenvektor mit Zeilenweisen Schlussfolgerungen an
    %S enthält jetzt Zeilenweise die Szenarien, inklusive der Conclusions in der letzten Spalte
    SMatrix = [SMatrix, SC];
    
	%Generiere Items (Jedes Item ist Bezeichung, z.B. (111|1), für ein Szenario) für Plot-Dropdownlist
    %A ist String-Szenarientabelle, Einträge sind Werte aus S im Datentyp String
    A = double(SMatrix);
    A = string(A);
    %gehe jedes Szenario zeilweise durch
    for i = 1:scens
        s = "";
        %Füge Prämisseneinzelwerte zu einem String zusammen. 
        for j = 1:premises
            s = strcat(s,A(i,j));
        end
        %Füge | als trenner zwischen prämissen und konklusion ein
        s = strcat(s,"|");
        %Hänge konklusion des szenarios an.
                
        s = strcat(s,A(i,(premises+1)));
        %s ist string der form "ppppp|c" mit belegungen aus
        %szenario
        string_SMatrix(i) = s;
    end
    string_SMatrix=[string_SMatrix,"Weighted Aggregation"];   
end
