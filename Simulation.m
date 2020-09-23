classdef Simulation < handle
%   Simulation class. Describes a simulation, simulates and stores
%   simulated data in an Simulation object
%			
%   Code from the Bachelor Thesis: 'Quantitative Epistemologische Demokratie:
%   Simulation und Vergleich der konklusions- und der Prämissenbasierten
%   Judgment Aggregation, Lukas Bosch, KIT
    
    properties
       %values from GUI
       name
       premises
       voters
       resultRelation
       runs
       pIterations
       %Correlation coefficient (All Covariances are equal)
       corr
       %premise prior probabilities (prob(1) = P(A=1), prob(2) = P(B=1) ...)
       priors
       %number of scenarios
       scens
       
       %Logical matrix, each row containis a szenario, last column is for conclusion
       SMatrix;
       %each line is a szenario representation in string form (e.g.(111|1) means all three premises and the conclusion are true (1)
       string_SMatrix
       %Scenario Structure-Array, each element is of type Scenario, 
       Scenarios;
       %Resultvector, sum of weighted scenarioresults, "Weighted Aggregation"
       Result
       wb
    end
    
    methods
        
        %%Konstructor
        function obj = Simulation(premises, voters, resultRelation, corr, runs, pIterations, name)

            "Version 4.0 "
            %set variables
            obj.name = name;
            obj.premises = premises;
            obj.scens = 2^premises;
            obj.voters = voters;
            obj.resultRelation = resultRelation;
            obj.corr = corr;
            obj.runs = runs;
            obj.pIterations = pIterations;
            %generate scenarios (SMatrix, obj.string_SMatrix)
            [obj.SMatrix, obj.string_SMatrix] = genScens(premises,resultRelation);
            obj.SMatrix
            %obj.string_SMatrix
            
            %Creating Scenario Structure
            for i=1:obj.scens
                obj.Scenarios(i).premises = premises;
                %1D Matrix with "Facts" (truth-values pf specific scenario);
                obj.Scenarios(i).Facts = obj.SMatrix(i,:);
                %PremiseFacts = Facts without last element (without conclusion). needed as Input-argument for signaltoVote
                obj.Scenarios(i).PFacts = obj.Scenarios(i).Facts(1:premises);
                %weighting factor w is calculated by setter method from gui
                %obj.Scenarios(i).w
                obj.Scenarios(i).goodCBM = 0;
                obj.Scenarios(i).goodPBM = 0;
                %resultvector
                obj.Scenarios(i).SResult = [];
            end
        end 
        
        %Sets prior probabilities. called from GUI
        function this = setPriors(this, premiseProbabilities)
            this.priors = premiseProbabilities;
        end
        
        %Sets scenario weights. called from GUI
        function this = setScenarioWeights(this)
           for i = 1:this.scens
               this.Scenarios(i).w = calcScenarioWeight(this.priors, this.Scenarios(i).Facts);
           end
        end
        
        %calculates weighted result for "Weighted Aggregation" Plot. called from GUI
        function this = calcWeightedResult(this)
           this.Result = zeros(length(this.pIterations),2);
           for i = 1:this.scens
              SR = this.Scenarios(i).SResult;
              this.Result = this.Result + (this.Scenarios(i).w * SR);
           end
        end
        
        %simulates this Simulation object for its parameters. called from
        %GUI
        function [this, finished] = simulate(this)
            
            %Progress Bar
            this.wb = waitbar(0, "", "Name", "Simulation Progress", 'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');%, "windowstyle", "modal");
            setappdata(this.wb, "canceling", 0);
            frames = java.awt.Frame.getFrames();
            frames(end).setAlwaysOnTop(1);
            
            %Declaration of symbolic variables
            syms A B C D E F G H I J K L M N O
            %generates symbolic expression (e.g. A&B) from resultRelation (e.g. "A&B")
            conclusion(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O) = str2sym(this.resultRelation);
            %generate anonymous function from symbolic expression
            conclude = matlabFunction(conclusion, "Vars",{[A B C D E F G H I J K L M N O]});
            
            %Pre-initialization of frequently used variables
            %sizeVC and sizeu are containing the dimensions of matrices
            %Voter Conclusion VC, mean u and majority maj
            sizeVC = [this.voters, 1];
            sizeu = [1, this.premises];
            VC = zeros(sizeVC);
            VC = logical(VC);
            maj = zeros(sizeu);
            u = zeros(sizeu);           

            length_pI = length(this.pIterations)
            for pIndex = 1:length_pI
                current_p = this.pIterations(pIndex);
                %check for cancel
                if getappdata(this.wb, "canceling")
                    finished = 0;
                    break
                end
                %update waitbar
                this.wb = waitbar(current_p,this.wb);
                
                %create variances and mean-vector from p
                mu(1:this.premises) = current_p;
                Var = mu.*(1-mu);
                %Calculate covarianze from correlation (Cov12 = corr12 * Var1 * Var2)
                %Variances are the same for all premises (Var(1) = Var(2) = ... = Var(n))
                %cova = this.corr * sqrt(Var(1)) * sqrt *(Var(2));
                cova = this.corr * Var(1);
                Cov = diag(Var);
                for i = 1:this.premises
                    for j = 1:this.premises
                        if i ~= j
                            Cov(i,j) = cova;
                        end
                    end
                end
                
                %Calculate Gauss u and Cov
                try
                    [Gu,GCov] = findLatentGaussian(mu,Cov);
                catch
                    delete(this.wb);
                    error('A joint Bernoulli distribution with the given covariance matrix does not exist!')
                end
                    
                %GCov = double(GCov);
                %Check for positive definiteness
                [chol_R, chol_p] = chol(GCov);
                if chol_p > 0
                    warning(['Covariance matrix of the latent Gaussian has at least one negative eigenvalue. ' ...
                                'Applying Higham-Correction (see help higham).'])
                    Cov
                    GCov = higham(GCov,1e-10,1e5)
                end
                
                %Monte-Carlo loop
                for k = 1:this.runs
                    
                    %Generate Signals (true or false / 1 or 0)
                    signals = sampleDichGauss01(Gu,GCov,this.voters, 1);
                    %Signals = transpose(Signals);
                    for i = 1:this.scens
                        %SizePfacts = size(this.Scenarios(i).PFacts)
                        %PFacts = this.Scenarios(i).PFacts
                        vote = abs(this.Scenarios(i).PFacts-(1-signals));
                        %calls majority. Maj ist vector containing the majority for each premise                   
                        if this.voters > 1
                            maj = majority(vote, this.voters);
                        else
                            maj = vote;
                        end
                        
                        %VC contains the conclusion of each voter
                        vc = conclude(vote);
                        %Conclusion-based majority/Premise-based majority
                        CBM = majority(vc, this.voters);
                        PBM = conclude(maj);
     
                        %Update goodCBM and goodPBM counter in Scenario object
                        if CBM == this.Scenarios(i).Facts(this.premises+1)                                                     
                            this.Scenarios(i).goodCBM = this.Scenarios(i).goodCBM + 1;
                        end
                        if PBM == this.Scenarios(i).Facts(this.premises+1)
                            this.Scenarios(i).goodPBM = this.Scenarios(i).goodPBM + 1;
                        end
                    end
                end                
                finished = 1;
                
                %End MC-loop. 
                %saving results for current_p for each Scenario in its SResult-Matrix. 
                %Reset of goodCBM and goodPBM counter because of new p in next p iteration
                for i=1:this.scens
                    %saving
                    Res = [this.Scenarios(i).goodCBM/this.runs, this.Scenarios(i).goodPBM/this.runs];
                    this.Scenarios(i).SResult = [this.Scenarios(i).SResult;Res];
                    %counter reset
                    this.Scenarios(i).goodCBM = 0;
                    this.Scenarios(i).goodPBM = 0;
                end
            end            
            delete(this.wb);
        end
    end
end


