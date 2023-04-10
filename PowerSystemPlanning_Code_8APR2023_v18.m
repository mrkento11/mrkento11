clear; clc;
format shortG
import java.text.*;     v = DecimalFormat;  % To use when printing results for easier readibility.
% Keep Track of time
initial = datetime;

% Directory
Database_name  = "Benchmark6bus";           %% Change name depending on which information we want to work with.
Route = "C:\Users\newsu\Documents\MATLAB" +"\\"+ Database_name;           %% The master directory will be the place where this code is stored.
cd(Route);                                  %% Change the directory to aim to the corresponding database.                     

% load demand into a variable
demand_base = readmatrix("Demand.xlsx");
opt_year    = "2050";                       %% opt_year must be a string.
fst_study_horizon = 2021;
last_study_horizon = str2double(opt_year);  %% This can definitely be more elegant.

%% Global information
crit_time = 0.01;                           %% duration of first demand block
blocks    = 20;     
interest  = 0.1;                            %% interest rate for projects. [0-1].

%% cluster demand
[base_LDC,it] = clustering_block(demand_base,crit_time,blocks);
disp("The amount of iterations for the cluster function to converge is: "+num2str(it))
disp("The Power of each cluster is: ")
% for readibility:
for i = 1:height(base_LDC)
    disp(char(v.format(base_LDC(i,2))) + " MW")     %% the representative power of each block
end
newline;
disp("The duration of each cluster is: ")
for i = 1:height(base_LDC)
    disp(char(v.format(base_LDC(i,1))) + " %")     %% the duration of each block
end
newline;

%% Create LDC Curve for base year 
Curve_LDC_base = [];
for t=1:blocks
    Curve_LDC_base = [Curve_LDC_base;ones(round(length(demand_base)*base_LDC(t,1)),1)*base_LDC(t,2)];
end

% Plot the LDC vs chronological demand for year 2019
graph_LDC_Smooth("2019",Curve_LDC_base,demand_base); % We use the function because we intend to use it several times
                                                     % OBSERVATION: Right now the figure is turned on. But just comment in case you don't want to see it.

%% Apply demand growth and obtain a set of new set of demands for the following years
% using information from the document: Generation Indicative Expansion Plan
% 2022-2031 Honduran Independent System Operator

% Invoke function to have tables with LDC and Smooth Demand
[Global_LDC_Demand,Global_Chronological_demand] = global_demands(Curve_LDC_base,demand_base);

% Graph for year to be optimized
% graph_LDC_Smooth(opt_year,Global_LDC_Demand.(opt_year),Global_Chronological_demand.(opt_year));

%% Get the hour-block map and the equivalent generation for all the blocks
% Hour-block map
hbm = hourblockmap(demand_base,base_LDC,opt_year,blocks);
% Equivalent generation for renewables for each block
[renew_gen_block, renewables_gen] = RENEW_GEN_BLK(hbm,blocks);

%% Create the Long Term Load and store it in a Table for optimization
% Obtain the LDC distribution for the last year of study
[opt_LDC,it] = clustering_block(Global_Chronological_demand.(opt_year),crit_time,blocks);
% Create the GEP_Load
GEP_Load = calc_GEP_Load(Global_LDC_Demand,opt_LDC,fst_study_horizon,last_study_horizon);

%% Renewable Commitment quota
Renew_commitment = calc_renew_commitment(opt_year);  %% SUGGESTION: Perhaps we just need milestone years according to plan

%% Iterations counter
iterations       = 0;
% Load general information from Database
DBB_info = call_BDD_info(Database_name + ".xlsx");
%% Flexibility and Capacity indicators
% Firm Capacity Reserve
sOR              = 0.10*ones(1,last_study_horizon-fst_study_horizon+1);    % Reserve Initial used for the initial problem
% Ramping
ramp_initial    = ramp_calc([],iterations,DBB_info);
ramp_capacity   = ramp_initial*ones(1,last_study_horizon-fst_study_horizon+1);  % ramp capacity for every year within the study horizon
% Define Renew_Cap
Renew_Cap = renew_cap_calc([],[],iterations,DBB_info,GEP_Load,[],[]);

% We make sure that the directory in which we will gonna graph is empty. To do so, in case information exists, we delete completely the directory. 
% In case it doesn't exist. We create the directory.
Reports_Directory = Route + "\\Reports";
if not(isfolder(Reports_Directory))
    mkdir(Reports_Directory);
elseif isfolder(Reports_Directory)
    rmdir(Reports_Directory,'s');
    mkdir(Reports_Directory);
end

% Information for peak_Yearly_Installed capacity graph
flag_peakLoad_InstalledCapacity = 0;    % Create a flag that will be used for the Peak Demand per year vs Peak Installed Capacity per Iteration
yearly_Installed_Capacity = [];         % Initialize the vector at zero
legend_description = [];                % Initialize the vector as empty

%% Create optimizer objects
yalmip('clear');
GEP_optimizer = GEP_optimizer_creator(GEP_Load,blocks,Renew_commitment,renew_gen_block,DBB_info,opt_year,fst_study_horizon,interest);    % Higher hierarchy GEP optimizer object
% Define the number of test days for our different criterias
qt_days = 7;    % We define how long will be the vector for days to be tested depending all of the criteria
UCP_optimizer = UCP_optimizer_creator(qt_days*24,DBB_info);         % Verifier UCP optimizer object

%% Main loop
while true                  % Infinite loop that will only be broken if we get to the finish line.
    iterations = iterations + 1;
    disp("This is iteration: " + num2str(iterations));
    % Display Time    
    final = datetime;
    disp(datetime);

    %% First Stage: Generation Expansion Problem
    disp("First Stage: Generation Expansion Problem")
    close all hidden;        % Attempt to free up memory

    % Evaluate the GEP Optimizer:
    [Investment,Output_G_exist_GEP,Output_G_cand_GEP,Cost,GEP_optimizer] = GEP_evaluate(GEP_optimizer,sOR,ramp_capacity,DBB_info,Renew_Cap);
    
    % Graph Findings
    graph_GEP(Output_G_exist_GEP,Output_G_cand_GEP,iterations,fst_study_horizon,"off",GEP_Load,DBB_info,sOR,ramp_capacity);   %% Graph a Report with the information of the recent findings.
    
    %% Second stage: Unit-Commitment Stage    
    for t = 1:last_study_horizon-fst_study_horizon+1
        % Signal start of for loop
        disp("Year currently under study: "+num2str(fst_study_horizon+t-1));            %% Flag to keep track of code advancements.
        
        % Create demand:
        % 1. Identify days for demand:
        test_days = ident_days_to_test(Global_Chronological_demand,num2str(fst_study_horizon+t-1),renewables_gen,Investment(:,t),length(demand_base),qt_days,DBB_info);  

        % 2. Use identified days to obtain the demand:
        [Week_Chrono_LNL,Week_Chrono_MNL,Week_Chrono_LRG,Week_Chrono_MRG,Week_Chrono_MDM,Week_Chrono_RAN] = Load_Chrono_to_test(num2str(fst_study_horizon+t-1),Global_Chronological_demand,test_days);

        %% Unit Commitment Problem
        % General Feasibility variable to control the loop between GEP and UC
        Feasibility = 0;
        Aggregate_Curtailment = 0;
        % 1. Test the UC Problem for Most Demand
        close all hidden;        % Attempt to free up memory        
        % Evaluate the UCP optimizer
        [Output_G_cand,Output_G_exist,MDM_marker,UCP_optimizer,Case_Curtailment] = UCP_evaluate(renewables_gen,test_days.most_Demand,num2str(fst_study_horizon+t-1),Investment(:,t),Week_Chrono_MDM,UCP_optimizer,"Most Demand",DBB_info);
        % Decide feasibility in case feasible, graph.
        Feasibility = Feasibility + MDM_marker;
        Aggregate_Curtailment = Aggregate_Curtailment + Case_Curtailment;
        if MDM_marker == 0
            graph_UCP(Output_G_exist,Output_G_cand,sOR(t),num2str(fst_study_horizon+t-1),"Most Demand","off",DBB_info,ramp_capacity(t));
        else
            disp("Most Demand is not feasible.")
        end
        
        % 2. Test the UC Problem for Most Net Load
        close all hidden;        % Attempt to free up memory       
        % Evaluate the UCP optimizer
        [Output_G_cand,Output_G_exist,MNL_marker,UCP_optimizer,Case_Curtailment] = UCP_evaluate(renewables_gen,test_days.most_NetLoad,num2str(fst_study_horizon+t-1),Investment(:,t),Week_Chrono_MNL,UCP_optimizer,"Most Net Load",DBB_info);
        % Decide feasibility in case feasible, graph.
        Feasibility = Feasibility + MNL_marker;
        Aggregate_Curtailment = Aggregate_Curtailment + Case_Curtailment;
        if MNL_marker == 0
            graph_UCP(Output_G_exist,Output_G_cand,sOR(t),num2str(fst_study_horizon+t-1),"Most Net Load","off",DBB_info,ramp_capacity(t));
        else
            disp("Most Net Load is not feasible.")
        end

        % 3. Test the UC Problem for Least Renewable Generation
        close all hidden;        % Attempt to free up memory
        % Evaluate the UCP optimizer
        [Output_G_cand,Output_G_exist,LRG_marker,UCP_optimizer,Case_Curtailment] = UCP_evaluate(renewables_gen,test_days.least_Renew_gen,num2str(fst_study_horizon+t-1),Investment(:,t),Week_Chrono_LRG,UCP_optimizer,"Least Renewable Generation",DBB_info);
        % Decide feasibility in case feasible, graph.
        Feasibility = Feasibility + LRG_marker;
        Aggregate_Curtailment = Aggregate_Curtailment + Case_Curtailment;
        if LRG_marker == 0
            graph_UCP(Output_G_exist,Output_G_cand,sOR(t),num2str(fst_study_horizon+t-1),"Least Renewable Generation","off",DBB_info,ramp_capacity(t));
        else
            disp("Least Renewable Generation is not feasible.")
        end
        
        % 4. Test the UC Problem for Least Net Load
        close all hidden;        % Attempt to free up memory        
        % Evaluate the UCP optimizer
        [Output_G_cand,Output_G_exist,LNL_marker,UCP_optimizer,Case_Curtailment] = UCP_evaluate(renewables_gen,test_days.least_NetLoad,num2str(fst_study_horizon+t-1),Investment(:,t),Week_Chrono_LNL,UCP_optimizer,"Least Net Load",DBB_info);
        % Decide feasibility in case feasible, graph.
        Feasibility = Feasibility + LNL_marker;
        Aggregate_Curtailment = Aggregate_Curtailment + Case_Curtailment;
        if LNL_marker == 0
            graph_UCP(Output_G_exist,Output_G_cand,sOR(t),num2str(fst_study_horizon+t-1),"Least Net Load","off",DBB_info,ramp_capacity(t));
        else
            disp("Least Net Load is not feasible.")
        end

        % 5. Test the UC Problem for Most Renewable Generation
        close all hidden;        % Attempt to free up memory        
        % Evaluate the UCP optimizer
        [Output_G_cand,Output_G_exist,MRG_marker,UCP_optimizer,Case_Curtailment] = UCP_evaluate(renewables_gen,test_days.most_Renew_gen,num2str(fst_study_horizon+t-1),Investment(:,t),Week_Chrono_MRG,UCP_optimizer,"Most Renewable Generation",DBB_info);
        % Decide feasibility in case feasible, graph.
        Feasibility = Feasibility + MRG_marker;
        Aggregate_Curtailment = Aggregate_Curtailment + Case_Curtailment;
        if MRG_marker == 0
            graph_UCP(Output_G_exist,Output_G_cand,sOR(t),num2str(fst_study_horizon+t-1),"Most Renewable Generation","off",DBB_info,ramp_capacity(t));
        else
            disp("Most Renewable Generation is not feasible.")
        end

        % 6. Test the UC Problem for a Random Week
        close all hidden;        % Attempt to free up memory        
        % Evaluate the UCP optimizer
        [Output_G_cand,Output_G_exist,RAN_marker,UCP_optimizer,Case_Curtailment] = UCP_evaluate(renewables_gen,test_days.Random,num2str(fst_study_horizon+t-1),Investment(:,t),Week_Chrono_RAN,UCP_optimizer,"Random days",DBB_info);
        % Decide feasibility in case feasible, graph.
        Feasibility = Feasibility + RAN_marker;
        Aggregate_Curtailment = Aggregate_Curtailment + Case_Curtailment;
        if RAN_marker == 0
            graph_UCP(Output_G_exist,Output_G_cand,sOR(t),num2str(fst_study_horizon+t-1),"Random","off",DBB_info,ramp_capacity(t));
        else
            disp("Random days is not feasible.")
        end

        % Check if the Problem was feasible
        if Feasibility ~= 0         % If one or more of the previous UCP are not feasible, we will increase the amount of reserve, to try to make them feasible within the next iterations. 
            % Verify if we have enough Firm Power
            if MDM_marker ~= 0 || MNL_marker ~= 0 || LRG_marker ~= 0
                disp("System Operating Reserve is incremented from "+ num2str(sOR(t)*100) +"% to "+ num2str((sOR(t)+.1)*100)  +"%."+newline);
                sOR(t) = sOR(t) + 0.1;        
            end
            % 
            if LNL_marker ~= 0 || MRG_marker ~= 0
                Renew_Cap = 0.85*renew_cap_calc(Investment(:,t),renew_gen_block,iterations,DBB_info,GEP_Load,Renew_Cap,t);
                disp("The amount of renewable energy has been constrained.")
            end
            
            % Additional information from the UCP: Amount of curtailment
            Cases
            if Aggregate_Curtailment/sum([MDM_marker;MNL_marker;LRG_marker;LNL_marker;MRG_marker;RAN_marker] == 0) >= 0.2
                ramp_actual = ramp_calc(Investment(:,t),iterations,DBB_info);    % We give the input about, according to the last investment decision we made, we need at least that.
                ramp_capacity(t) = ramp_actual*1.5;                 % But we increase it by 50%.
                disp("The system ramping capability is increased from " + char(v.format(ramp_actual))+ " MW to "+ char(v.format(ramp_actual*1.5))+ " MW." +newline);
                final = datetime;
            end

            % Execute the graph_Peak_Invest function to keep track of the yearly Installed capacity
            [yearly_Installed_Capacity,legend_description] = graph_Peak_Invest(Investment,GEP_Load,flag_peakLoad_InstalledCapacity,yearly_Installed_Capacity,legend_description,renew_gen_block,blocks,DBB_info,fst_study_horizon+t-1,sOR(t),ramp_capacity(t));
            
            % We break the for loop
            final = datetime;
            break
        end
        % Additional information from the UCP: Amount of curtailment
        if Aggregate_Curtailment/sum([MDM_marker;MNL_marker;LRG_marker;LNL_marker;MRG_marker;RAN_marker] == 0) >= 0.2
            ramp_actual = ramp_calc(Investment(:,t),iterations,DBB_info);    % We give the input about, according to the last investment decision we made, we need at least that.
            ramp_capacity(t) = ramp_actual*1.5;                 % But we increase it by 50%.
            disp("The system ramping capability is increased from " + char(v.format(ramp_actual))+ " MW to "+ char(v.format(ramp_actual*1.5))+ " MW." +newline);
            final = datetime;

            % Execute the graph_Peak_Invest function to keep track of the yearly Installed capacity
            [yearly_Installed_Capacity,legend_description] = graph_Peak_Invest(Investment,GEP_Load,flag_peakLoad_InstalledCapacity,yearly_Installed_Capacity,legend_description,renew_gen_block,blocks,DBB_info,fst_study_horizon+t-1,sOR(t),ramp_capacity(t));
            break
        end

        disp("Year "+ num2str(fst_study_horizon+t-1) + " finished." +newline);            %% Flag to keep track of code advancements.
    end
    
    %% We get out from the while loop
    if and(t == last_study_horizon-fst_study_horizon+1, Feasibility == 0)
        flag_peakLoad_InstalledCapacity = flag_peakLoad_InstalledCapacity+1;
        [yearly_Installed_Capacity,legend_description] = graph_Peak_Invest(Investment,GEP_Load,flag_peakLoad_InstalledCapacity,yearly_Installed_Capacity,legend_description,renew_gen_block,blocks,DBB_info,fst_study_horizon+t-1,sOR(t),ramp_capacity(t));
        break
    end
end

final = datetime;
% Save the data to a mat file for review even after MatLab interface has been closed.
clear DBB_info;         % This is an nonessential variable, used only to handle in an object all the database information. However, this are know parameters stored initially in Excel files, so it is useless to store them in our final mat file. Other variables are nonessential, and could as well be deleted.
save(Reports_Directory +"\\"+ Database_name + ".mat");

%% Functions 

function [LDC_k,it] = clustering_block(demand_base,crit_time,blocks)
        % We will base most of the reasoning in here in the file: Pruebas Herramienta para determinar BH_ 12 MAR 2019
    
    % Info for clusters Block 1 = % of time for critical time Rest of the blocks:
    init_duraci = (1-crit_time)/(blocks-1);
    % Sort from greatest to smallets
    sorted_demand = sort(demand_base,'descend');
    % Look for centroids
    % Centroids and its Length
    Centroids   = zeros(blocks,1);
    Lengths     = zeros(blocks,1);
    % Initialization
    Lengths(1)   = round(length(demand_base)*crit_time); % Length or duration of block 1
    Centroids(1) = mean(sorted_demand(1:Lengths(1))); % This centroid represents the average of first duration length
    
    % Centroids and Duration of the rest is:
    for i=2:blocks
        if i == blocks
            Lengths(i)   = Lengths(i-1) + round(length(demand_base)*init_duraci)-1;     % We are accumulating the reference of the duration length
            Centroids(i) = mean(sorted_demand(Lengths(i-1)+1:Lengths(i)));              % we average those said duration lengths
        else
            Lengths(i)   = Lengths(i-1) + round(length(demand_base)*init_duraci);
            Centroids(i) = mean(sorted_demand(Lengths(i-1)+1:Lengths(i)));
        end
    end

    % Make a table of demand and distances to all Centroids
    Table = [];
    % Distances to each centroids
    for i =1:blocks
        Table(:,i) = abs(sorted_demand-Centroids(i));
    end

    % Assign flag to assign Centroid
    [~,Duration] = min(Table,[],2);
    % Count the duration of each Block
    [duraci,~] = groupcounts(Duration);
    
    % Update Lengths and Centroids
    Lengths(1) = duraci(1);
    Centroids(1) = mean(sorted_demand(1:Lengths(1)));
    % for loop to update blocks 2 till the finish
    for i = 2:blocks
        Lengths(i)   = Lengths(i-1) + duraci(i);
        Centroids(i) = mean(sorted_demand(Lengths(i-1)+1:Lengths(i)));
    end
    
    % Create Load Duration Curve information output
    Power = zeros(blocks,1);
    % peak demand will be preserved
    Power(1) = max(sorted_demand);
    % All the rest of the blocks = centroids
    for i=2:blocks
        Power(i) = Centroids(i);
    end


    % Correct First duration
    peak_excess = round((duraci(1) - round(length(demand_base)*crit_time))/(blocks-1));      % variable that will determine if we have an excess of number hours to adjust at peak demand
    duraci(1) = round(length(demand_base)*crit_time); % correct it back to the fixed duration block
    duraci(2:end) = duraci(2:end) + peak_excess;          % we will adjust with the duration of the rest of the blocks equally

    % Correct last Power
    Power(blocks) = abs((sum(sorted_demand)-Power(1:blocks-1)'*duraci(1:blocks-1)))/duraci(blocks);
    LDC = [duraci./8760,Power];
    

    % Start Counter for iterations
    it = 0;
    
    % Assign LDC_k to our current iteration
    LDC_k = LDC;
    LDC_k_1 = zeros(size(LDC_k,1),size(LDC_k,2));
    % While Loop to make convergence of 
    while sum(abs(LDC_k_1(:,2)-LDC_k(:,2))) >= blocks % the sum of all differences must be less than 1/20 MW times blocks
        LDC_k = LDC_k_1;
        % Update table of demand and distances to all Centroids
        for i =1:blocks
            Table(:,i) = abs(sorted_demand-Centroids(i));
        end
        
        % Assign flag to assign Centroid
        [~,Duration] = min(Table,[],2);
        % Count the duration of each Block
        [duraci,~] = groupcounts(Duration);

        % Update Lengths and Centroids
        Lengths(1) = duraci(1);
        Centroids(1) = mean(sorted_demand(1:Lengths(1)));
        % for loop to update blocks 2 till the finish
        for i = 2:blocks
            Lengths(i)   = Lengths(i-1) + duraci(i);
            Centroids(i) = mean(sorted_demand(Lengths(i-1)+1:Lengths(i)));
        end

        % Create Load Duration Curve information output
        Power = zeros(blocks,1);
        % peak demand will be preserved
        Power(1) = max(sorted_demand);
        % All the rest of the blocks = centroids except the last
        for i=2:blocks-1
            Power(i) = Centroids(i);
        end
    
        % Create Load Duration Curve information output
        Power = zeros(blocks,1);
        % peak demand will be preserved
        Power(1) = max(sorted_demand);
        % All the rest of the blocks = centroids except the last
        for i=2:blocks
            Power(i) = Centroids(i);
        end

        % Correct First duration
        peak_excess = round((duraci(1) - round(length(demand_base)*crit_time))/(blocks-1));      % variable that will determine if we have an excess of number hours to adjust at peak demand
        duraci(1) = round(length(demand_base)*crit_time); % correct it back to the fixed duration block
        duraci(2:end) = duraci(2:end) + peak_excess;          % we will adjust with the duration of the rest of the blocks equally

        % Correct last Power
        Power(blocks) = abs(sum(sorted_demand)-Power(1:blocks-1)'*duraci(1:blocks-1))/duraci(blocks);
        LDC_k_1 = [duraci./8760,Power];

        if it >= 100
            break
        else
            it = it +1;
        end
    end
    LDC_k = LDC_k_1;
end

function [Global_LDC_Info,Global_Chronological_Info] = global_demands(Curve_LDC,demand_base)
    % load growth into a variable
    demand_growth = readmatrix("DemandGrowth.csv");

    % Create Table with Info
    Global_LDC_Info = table(Curve_LDC);
    for ii = 1:length(demand_growth)
        % For the first year that is not 2019, we multiply against the first
        % growth factor
        if ii == 1
            current_demand = Curve_LDC*demand_growth(ii,2);
        % for the rest of the years, we multiply against the previous years we
        % just calculated
        else
            current_demand = current_demand*demand_growth(ii,2);
        end
        % We concatenate this into our Table
        Global_LDC_Info.(ii+1) = current_demand;
    end
    
    % We now will change the headers in our table
    % Initialize the table headers with year 2019
    headers = ["2019"];
    % For the rest we extract from the variable demand growth first column
    headers = [headers;num2str(demand_growth(:,1))];
    
    Global_LDC_Info.Properties.VariableNames = headers;
    
    % Now we repeat the procedure but for the chronological curve
    Global_Chronological_Info = table(demand_base);
    for ii = 1:length(demand_growth)
        % For the first year that is not 2019, we multiply against the first
        % growth factor
        if ii == 1
            current_demand = demand_base*demand_growth(ii,2);
        % for the rest of the years, we multiply against the previous years we
        % just calculated
        else
            current_demand = current_demand*demand_growth(ii,2);
        end
        % We concatenate this into our Table
        Global_Chronological_Info.(ii+1) = current_demand;
    end
    
    % We apply the same headers to the chronological demand
    Global_Chronological_Info.Properties.VariableNames = headers;
end

function graph_LDC_Smooth(year,Curve_LDC,chronological_demand)
    figure("Name","LDC vs Smooth Demand " + year,'NumberTitle','off','Visible',"on");
    
    hold on
    
    plot(Curve_LDC);
    plot(sort(chronological_demand,"descend"));
    ylim([0 max(chronological_demand)+200]);
    
    hold off
    legend('Load Duration Curve','Smooth Demand',"Location","best")
end

function C = hourblockmap(demand_base,base_LDC,opt_year,blocks)
    duraci = length(demand_base)*base_LDC(:,1); % duration of each block independently

    % Create a vector for duration aggregated. This means that that specific block goes from the previous number up to the number stored for that block. Ex: 1-> 613   2-> 3986. So block 2 goes from 614-3986
    duraci_agg = zeros(blocks,2);
    duraci_agg(1,1) = duraci(1);
    for i=2:blocks
        duraci_agg(i,1) = duraci_agg(i-1,1) + duraci(i,1);
    end
    % Now we will paste a second column with the corresponding block flag to it
    blk_id = 1;
    for i = 1:length(duraci_agg)
        duraci_agg(i,2) = blk_id;
        blk_id = blk_id + 1;
    end
    
    % Long vector for flag about hourblock map  
    blk_flag = zeros(length(demand_base),1); % this vector has the same length as the yearly demand
    blk_id   = 1;                            % we use again a flag to mark the varible blk_flag
    for ii = 1:length(demand_base)           % The logic is that the for goes through the whole blk_flag, if that row number is equal or less than the duraci_agg() row number, then it prints the second column
        if ii <= duraci_agg(blk_id,1)
            blk_flag(ii) = duraci_agg(blk_id,2);
            if ii == duraci_agg(blk_id,1)    % If it reaches the extreme point, we go to the next row of duraci_agg by incrementing blk_id + 1.
                blk_id = blk_id + 1;
            end
        end
    end
    
    % Use blocks flag in ordered demand base
    year = str2double(opt_year);
    % First make a first date
    t = datetime([year,1,1]);
    % Second, converti it to a number
    t = datenum(t);
    % Third, assign to time_stamps(1) as initial date
    time_stamps = NaN(length(demand_base),1);
    time_stamps(1) = t;
    % fourth, fill up the correct numbers for all the long vector time
    for ii = 2:length(time_stamps)
        time_stamps(ii) = time_stamps(ii-1) + 1/24; % for the length in time_stamps we add one hour
    end
    
    % Finally, we convert back to date format, while at the same time round
    % to the nearest hour
    time_stamps = dateshift(datetime(time_stamps,'ConvertFrom','datenum'),'start','hour','nearest');
    
    
    % Create a table with the date and the base demand
    A = table(time_stamps,demand_base);
    % Sort the table upon the demand
    A = sortrows(A,'demand_base','descend');
    % take that sorted table and the time stampos to attach it to the block
    % flags
    B = table(A.time_stamps,blk_flag);
    % Just rename variables
    B.Properties.VariableNames = ["time_stamps","block_flag"];
    % Sort again by date in chronological manner
    C = sortrows(B,"time_stamps","ascend");
end

function [renew_gen_block, renewables_gen] = RENEW_GEN_BLK(hbm,blocks)
    % This function takes in the hour-block map (hbm) and %blocks (blocks), and delivers the block renewable output generation in block
    % (renew_gen_block) and the hourly renewable output generation in matrix form (renewables_gen).
    % The function is divided into three steps: 
    % 1. Read the hourly generation from "RenewablesScenarios.xlsb" taken from Renewables Ninja
    % 2. Add hourblock flag to the table to make the average if regarding this criteria.
    % 3. Make an average if depending on the hour block map. Documentation: https://www.mathworks.com/matlabcentral/answers/405518-find-the-mean-of-values-in-a-column-when-a-condition-is-met

    % 1. Take data as input
    renewables_gen = readtable("RenewableScenario.xlsx",'ReadVariableNames',true,"UseExcel",false);      % The file has the variable names in the first line, so it is needed the ReadVariablesNames == true.
    RenewableScenariosName = renewables_gen.Properties.VariableNames(2:end);            % not taking into account the first one because that one is the time_stamp. Just considering the RenewableScenarios.

    % 2. Add hourblock flag to the table to make the average if regarding this criteria.
    % Add hourblock flags to the table
    % hbm_flag = table2array(hbm(:,2));
    renewables_gen = addvars(renewables_gen,hbm.block_flag,'After','time_stamp','NewVariableNames',['block_flag']);
    
    % 3. Make average if for all renewables
    renew_gen_block = zeros(blocks,width(renewables_gen)-2);                % Matrix dimension is {blocks, renewable sceanarios (width_renewables_gen - column for hbm and column for time_stamp)}
    renewables_gen_cpy = removevars(renewables_gen,"time_stamp");           
    renewables_gen_cpy = table2array(renewables_gen_cpy);
    for ii = 2:width(renewables_gen_cpy)
        for jj = 1:blocks   
            renew_gen_block(jj,ii-1) = mean(renewables_gen_cpy(renewables_gen_cpy(:,1) == jj, ii));   % Mean Of Column ii For Column 'hbm_flag' = jj
        end                                                                                           % jj runs from 1:5 as blocks
    end                                                                                               % ii runs from third to last column un renewables_gen
    
    % Convert renewewables generation block to table
    renew_gen_block = array2table(renew_gen_block);
    renew_gen_block.Properties.VariableNames = [RenewableScenariosName'];
    
    
    % C_mean = mean(M(M(:,1) == someValue, 3))                % Mean Of Column #3 For Column #1 = %somevalue 
    % This last string is stored to preserve an example on how to possible
    % perform an averageif on a table. Still haven't figured it out.
    % Another possible example is:
    % renew_gen_block = mean(renewables_gen{:,vars}) this line gives an
    % absolute average, but not quite sure how to make the average if for all
    % vars.
end

function GEP_optimizer = GEP_optimizer_creator(GEP_Load,blocks,Renew_commitment,renew_gen_block,DBB_info,opt_year,fst_study_horizon,interest)
    %% Call Input Information about existing and candidate generators
    Thermal_exist_data = DBB_info.Thermal_exist_data;
    EcoThermal_exist_data = DBB_info.EcoThermal_exist_data;
    Renewables_exist_data = DBB_info.Renewables_exist_data;    
    Thermal_cand_data = DBB_info.Thermal_cand_data;
    EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
    Renewables_cand_data = DBB_info.Renewables_cand_data;
    ESS_cand_data = DBB_info.ESS_cand_data;
    Fuel_price_data = DBB_info.Fuel_price_data;
    Flags = DBB_info.Flags;

    %% General Information for the whole mathematical problem. This data shouldn't change across the years.
    % Note: The renewable information might change in the future. As of now, the problem will be deterministic and the generation of
    % renewable will be the same across all years. However, to make it stochastic, we will need to create a function that chooses how
    % renewables will generate energy across each year.

    % Existing Generator Bounds
    Exist_Upper = [Thermal_exist_data.UpperLimit; EcoThermal_exist_data.UpperLimit; Renewables_exist_data.UpperLimit]; % Upper Generation Bounds Existing
    
    % Candidate Generator Bounds
    Cand_Upper = [Thermal_cand_data.UpperLimit; EcoThermal_cand_data.UpperLimit; Renewables_cand_data.UpperLimit; ESS_cand_data.UpperLimit]; % Upper Generation Bounds Candidate
    
    % Max Installation Bounds: Given that decision variable is integer, we set
    % a limit on the total install capacity coming from each source
    MaxInstall = [Thermal_cand_data.MaxInstall; EcoThermal_cand_data.MaxInstall; Renewables_cand_data.MaxInstall; ESS_cand_data.MaxInstall];
    
    % Load
    PLoad = GEP_Load{:,2:end}; % Power Instances will be equal to the representing blocks
    
    % Renewable generation quota
    % Peaking Information: According to information from paper couple long-term
    % and short-term
    vCap_exist= Exist_Upper;    % the Capacity of Installed plants is the same as the Upper Limit of Existing plants
    % Peak Factor coefficient factor contribution for existing generators
    % Renewables Peak Factor coefficient factor
    Exist_renewables_capacity_factor = ones(blocks,height(Renewables_exist_data));
    for i = 1:height(Renewables_exist_data)
        Exist_renewables_capacity_factor(:,i) = [renew_gen_block.(Renewables_exist_data.RenewableScenario{i})];
    end
    % Take out the capacity factor of renewables for peak demand
    pPF_exist = [Thermal_exist_data.HistoricalAvailability;EcoThermal_exist_data.HistoricalAvailability;Exist_renewables_capacity_factor(1,:)']; % Peak factor of all technologies in peak demand.
    
    % Peak Factor coefficient factor contribution for candidate generators
    Cand_renewables_capacity_factor = ones(blocks,height(Renewables_cand_data));
    for i = 1:height(Renewables_cand_data)
        Cand_renewables_capacity_factor(:,i) = [renew_gen_block.(string(Renewables_cand_data.RenewableScenario{i}))];
    end
    % Take out the capacity factor of renewables for peak demand
    pPF_cand = [Thermal_cand_data.HistoricalAvailability;EcoThermal_cand_data.HistoricalAvailability;Cand_renewables_capacity_factor(1,:)';ones(height(ESS_cand_data),1)]; % Peak factor of all technologies in peak demand.        We define the factor of Batteries as able to deliver all of their installed capacity.
    
    %% Mathematical Formulation
    % Decision variables
    xit            = intvar(Flags.quant_cand,width(GEP_Load)-1,'full');                     % Integer variables to define the installed capacity needed. Variable is 2-D: [quantity of type of technology decision, years]. The decision variables doesn't change per block.
    yit            = intvar(Flags.quant_cand,width(GEP_Load)-1,'full');                     % Integer variable used to measure Investment. Not a decision variable. Dependent variable.
    wit            = intvar(Flags.quant_cand,width(GEP_Load)-1,'full');                     % Integer variables used to measure Amortization variables 
    Exist_git      = sdpvar(Flags.quant_exist,blocks,width(GEP_Load)-1,'full');             % continuous variable to define generation output from exist generation. Variable is 3-D, for same reason as Cand_git.
    Cand_git       = sdpvar(Flags.flag_renewables_cand,blocks,width(GEP_Load)-1,'full');    % continuous variable to define generation output from future generation. Variables is 3-D: [quantity of candidate generators, block, year]. There should be an independent value for each gen i, block b, and year t.
    q_Charge_itb   = sdpvar(height(ESS_cand_data),blocks,width(GEP_Load)-1,'full');         % continuous variable to define discharge from candidate batteries.
    q_Discharge_itb = sdpvar(height(ESS_cand_data),blocks,width(GEP_Load)-1,'full');        % continuous variable to define discharge from candidate batteries.
    % Decision variables for curtailment
    Curt_exist_git = sdpvar(Flags.flag_renewables_exist-Flags.flag_EcoThermal_exist,blocks,width(GEP_Load)-1,'full');             % continuous variable to define curtailment from renewables. Variable is 3-D, for same reason as Exist_git.
    Curt_cand_git = sdpvar(Flags.flag_renewables_cand-Flags.flag_EcoThermal_cand,blocks,width(GEP_Load)-1,'full');             % continuous variable to define curtailment from renewables. Variable is 3-D, for same reason as Exist_git.

    % Nomenclature:
    % i: Generators. Will be used both for aggregated generators in Existing or Candidate decision variable object.
    % b: blocks.
    % t: years.
    
    % Parameters: parameters that make certain constraints change over time.
    sOR             = sdpvar(1,width(GEP_Load)-1,'full');       % System Operating Reserve parameter
    ramp_capacity   = sdpvar(1,width(GEP_Load)-1,'full');       % ramping capacity for each year.
    Renew_Cap       = sdpvar(1,width(GEP_Load)-1,'full');       % Renewable Energy Cap for each year to prevent too much curtailment.

    %% Constraints
    % General
    Constraints = [];
    % Generator Bounds: Min and Max Constraints
    % Existing Generator Bounds
    for t = 1:width(GEP_Load)-1
        for b = 1:blocks
            Constraints = Constraints + [(0 <= Exist_git(:,b,t) <= Exist_Upper) : 'Existing Generator Bounds'];       % given that we are not considering on and off conditions for the machines, the lower bound is zero.
        end
    end
    % Candidate Generator Bounds: All technologies including batteries
    for t = 1:width(GEP_Load)-1
        for b = 1:blocks                % vector of zeros size as long as Git_cand + Batteries <= [Git_cand;Batteries]              <= decisionvariable*UpperLimit
            Constraints = Constraints + [(0 <= [Cand_git(:,b,t);q_Discharge_itb(:,b,t)] <= xit(:,t).*Cand_Upper): 'Candidate Generator Bounds']; % given that we are not considering on and off conditions for the machines, the lower bound is zero.
            % Given that iterates over the same loops, we will add the battery charge constraint in here
            Constraints = Constraints + [(0 <= q_Charge_itb(:,b,t) <= xit(Flags.flag_renewables_cand+1:end,t).*ESS_cand_data.UpperLimit): 'Battery Bounds'];       % charge constraint
        end                            % vector of zeros size as big as Batteries <=  Batteries charge <= decision variable*UpperLimit
    end

    % Power Balance Constraint
    for t = 1:width(PLoad)
        for b = 1:blocks
            Constraints = Constraints + [(sum(Exist_git(:,b,t)) + sum(Cand_git(:,b,t)) + sum(q_Discharge_itb(:,b,t)) == PLoad(b,t) + sum(q_Charge_itb(:,b,t))): 'Power Balance Constraint'];          % Fix batteries in here
        end     % This constraint reads: Sum(Exist_Git) + Sum(Candi_Git)              + batteries discharge         = Load for that specific block + batteries charge
    end
    
    % Max Installation constraint: perhaps we cannot build more than a certain amount of technology given area constraints
    for t = 1:width(GEP_Load)-1
        Constraints = Constraints + [[xit(:,t).*Cand_Upper <= MaxInstall] : 'Max Install'];
    end

    % Peaking Equation: Firm Power Security
    for t = 1:width(GEP_Load)-1
        Constraints = Constraints + [[(pPF_exist'*vCap_exist + (xit(:,t).*pPF_cand)'*Cand_Upper) >= PLoad(1,t)*(1+sOR(t))]: 'Peaking Equation'];          %sOR: system Operating Reserve
    end                             % Reserve * (pPF_exist*InstalledCap + decision variable*pPF_cand*MaxCapacity) >= Peak demand at each year.

    % Ramping Capacity:
    for t = 1:width(GEP_Load)-1     % Ramp capacity at each year must be at least the RampUp of Existing Thermal, EcoThermal + RampUp of Candidate Thermal,EcoThermal + Installed Capacity of Batteries. We consider the installed capacity of Batteries because they can act immediately.
        Constraints = Constraints + [(ramp_capacity(t) <= sum(Thermal_exist_data.RampUp) + sum(EcoThermal_exist_data.RampUp) + Thermal_cand_data.RampUp'*xit(1:Flags.flag_thermal_cand,t) + EcoThermal_cand_data.RampUp'*xit(Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand,t) + ESS_cand_data.MaxInstall'*xit(Flags.flag_renewables_cand+1:Flags.flag_ESS_cand,t)): 'Ramping constraint'];
    end

    % Cap on the amount of energy that can come from Renewables
    for t = 1:width(GEP_Load)-1     % sum(Exist_Upper*AvailableResource) + sum(xit*Exist_Upper*AvailableResource) <= Renew_Cap
        LHS = 0;
        % Existing Renewables
        counter = 1;        
        for i = Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist
            LHS = LHS + sum(Exist_Upper(i).*renew_gen_block.(string(Renewables_exist_data.RenewableScenario{counter})));
            counter = counter + 1;
        end
        % Candidate Renewables
        counter = 1;        
        for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
            LHS = LHS + sum(xit(i,t)*Cand_Upper(i).*renew_gen_block.(string(Renewables_cand_data.RenewableScenario{counter})));
            counter = counter + 1;
        end
        % Actual Constraint
        Constraints = Constraints + [LHS <= Renew_Cap(t)];
    end
    
    %% Renewables generation resource: This constraint is read that the generation from renwables needs to be at most the available resource
    % Existing renewables
    for t = 1:width(GEP_Load)-1    
        counter = 1; % fictitious variable to start the count for the length of Solar
        for i = Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist   %% i starts where the thermal generators ends, and for an amount of solar generators as defined by the file Exist_Solar.
            Constraints = Constraints + [(Exist_git(i,:,t)' <= Exist_Upper(i).*renew_gen_block.(string(Renewables_exist_data.RenewableScenario{counter}))): 'Existing Renewable resource'];
            counter = counter + 1;      % This constraint reads: Git <= InstalledCapacity*resource for that time and block.
        end
    end
    
    % Candidate renewables
    for t = 1:width(GEP_Load)-1
        counter = 1; % fictitious variable to start the count for the length of Solar
        for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand     %% i starts where thermal generators end, and for an amount of solar as defined in Cand_Solar.
            Constraints = Constraints + [(Cand_git(i,:,t)' <= xit(i,t)*Cand_Upper(i).*renew_gen_block.(string(Renewables_cand_data.RenewableScenario{counter}))): 'Existing Renewable resource'];
            counter = counter + 1;      % This constraint reads: Git <= InstalledCapacity*resource for that time and block.
        end
    end

    % Add constraint to make solar and wind to not grow in an uneven way
    % Find where the data for wind starts.
    first_wind_exist = find(Renewables_exist_data.Tech == "Wind",1);
    first_wind_cand = find(Renewables_cand_data.Tech == "Wind",1);  
    for t = 1:width(GEP_Load)-1
        Constraints = Constraints + [0 <= sum(Exist_git(Flags.flag_EcoThermal_exist+1:Flags.flag_EcoThermal_exist+first_wind_exist-1,:,t)*(8760*GEP_Load.LD)) + sum(Cand_git(Flags.flag_EcoThermal_cand+1:Flags.flag_EcoThermal_cand+first_wind_cand-1,:,t)*(8760*GEP_Load.LD)) <= 3*(sum(Exist_git(Flags.flag_EcoThermal_exist+first_wind_exist:Flags.flag_renewables_exist,:,t)*(8760*GEP_Load.LD)) + sum(Cand_git(Flags.flag_EcoThermal_cand+first_wind_cand:Flags.flag_renewables_cand,:,t)*(8760*GEP_Load.LD)))];
        Constraints = Constraints + [0 <= sum(Exist_git(Flags.flag_EcoThermal_exist+first_wind_exist:Flags.flag_renewables_exist,:,t)*(8760*GEP_Load.LD)) + sum(Cand_git(Flags.flag_EcoThermal_cand+first_wind_cand:Flags.flag_renewables_cand,:,t)*(8760*GEP_Load.LD)) <= 3*(sum(Exist_git(Flags.flag_EcoThermal_exist+1:Flags.flag_EcoThermal_exist+first_wind_exist-1,:,t)*(8760*GEP_Load.LD)) + sum(Cand_git(Flags.flag_EcoThermal_cand+1:Flags.flag_EcoThermal_cand+first_wind_cand-1,:,t)*(8760*GEP_Load.LD)))];
    end
    
    % Renewable quota contraint. It reads: Sum(Cand_Git(StartAtEcoThermal:end,:,final year)*                          *(vector of LoadDuration)) + Sum(Exist_Git(StartAtEcoThermak:end,:,final year)*                                 *(vector of LoadDuration)) >= Renewable commitment*8760*LoadDuration*Clusters;
    Constraints = Constraints + [[sum(Cand_git(Flags.flag_thermal_cand+1:Flags.flag_renewables_cand,:,str2double(opt_year)-fst_study_horizon+1)*(8760*GEP_Load.LD)) + sum(Exist_git(Flags.flag_thermal_exist+1:Flags.flag_renewables_exist,:,str2double(opt_year)-fst_study_horizon+1)*(8760*GEP_Load.LD)) >= Renew_commitment*8760*GEP_Load.LD'*GEP_Load.(opt_year)]: 'Renewable Quota']; 

    % Measure Renewable Curtailment
    % Existing renewables
    for t = 1:width(GEP_Load)-1    
        counter = 1;                % We create this counter because the index i is different when pointing in the Renewable Scenario list, than when pointing towards the position inside the Exist_Upper and Exist_git matrices.
        for i = Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist
            Constraints = Constraints + [[Curt_exist_git(counter,:,t)' == Exist_Upper(i).*renew_gen_block.(string(Renewables_exist_data.RenewableScenario{counter})) - Exist_git(i,:,t)']: 'Existing Renewable Curtailment'];
            counter = counter + 1;      % Constraint reads: Curtailment == Renewable Energy Available - Actual Renewable Energy Used
        end
    end

    for t = 1:width(GEP_Load)-1
        counter = 1;                % We create this counter because the index i is different when pointing in the Renewable Scenario list, than when pointing towards the position inside the Cand_Upper and Cand_git matrices.
        for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
            Constraints = Constraints + [[Curt_cand_git(counter,:,t)' == xit(i,t)*Cand_Upper(i).*renew_gen_block.(string(Renewables_cand_data.RenewableScenario{counter})) - Cand_git(i,:,t)']: 'Candidate Renewable Curtailment'];
            counter = counter + 1;      % Constraint reads: Curtailment == Renewable Energy Available - Actual Renewable Energy Used
        end
    end

    % Maximum Curtailment per year
    % Existing Generators
    for t = 1:width(GEP_Load)-1
        counter = 1;
        for i = Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist
            % This constraint reads:    sum(Curtailment of Renewable i at year t) <= 0.05*Installed Capacity*Actual resource. So, we are saying that we can only curtail 5% of the resource at Planning stage.
            Constraints = Constraints + [sum(Curt_exist_git(counter,:,t)) <= 0.03*sum(Exist_Upper(i)*renew_gen_block.(string(Renewables_exist_data.RenewableScenario{counter})))];
            counter = counter + 1;
        end
    end
    % Candidate Generators
    for t = 1:width(GEP_Load)-1
        counter = 1;
        for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
            Constraints = Constraints + [sum(Curt_cand_git(counter,:,t)) <= 0.03*sum(xit(i,t)*Cand_git(i)*renew_gen_block.(string(Renewables_cand_data.RenewableScenario{counter})))]; 
        end
    end


    %% Minimum year for investments: All Candidates
    % Thermal
    for t = 1:width(GEP_Load)-1
        for i = 1:Flags.flag_thermal_cand
            if Thermal_cand_data.MinYear(i) > fst_study_horizon +t -1
                Constraints = Constraints + [(xit(i,t) == 0): 'Thermal minimum year constraint'];
            end
        end
    end
    % EcoThermal
    for t = 1:width(GEP_Load)-1
        counter = 1;
        for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
            if Thermal_cand_data.MinYear(counter) > fst_study_horizon +t -1
                Constraints = Constraints + [(xit(i,t) == 0): 'EcoThermal minimum year constraint'];
            end
            counter = counter +1;
        end
    end
    % Renewables
    for t = 1:width(GEP_Load)-1
        counter = 1;
        for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
            if Thermal_cand_data.MinYear(counter) > fst_study_horizon +t -1
                Constraints = Constraints + [(xit(i,t) == 0): 'Renewables minimum year constraint'];
            end
            counter = counter +1;
        end
    end
    % Energy Storage Systems
    for t = 1:width(GEP_Load)-1
        counter = 1;
        for i = Flags.flag_renewables_cand+1:Flags.flag_ESS_cand
            if Thermal_cand_data.MinYear(counter) > fst_study_horizon +t -1
                Constraints = Constraints + [(xit(i,t) == 0): 'ESS minimum year constraint'];
            end
            counter = counter +1;
        end
    end
    
    %% Battery Additional Constraints
    % Battery Energy Conservation Constraint
    for t = 1:width(GEP_Load)-1
        for i = 1:height(ESS_cand_data) % the charged and discharged energy decision variables do not need a counter flag because it only is defined by the amount of batteries in the system
            Constraints = Constraints + [(q_Discharge_itb(i,:,t)*GEP_Load.LD*8760 == ESS_cand_data.Efficiency(i)*8760*(q_Charge_itb(i,:,t)*GEP_Load.LD)): 'Batteries Energy conservation constraint'];
        end
    end
    
    % Maximum energy discharged per year
    for t = 1:width(GEP_Load)-1
        for i = 1:height(ESS_cand_data)
            Constraints = Constraints + [(q_Discharge_itb(i,:,t)*GEP_Load.LD*8760 <= 365*ESS_cand_data.Duration(i)*ESS_cand_data.UpperLimit(i)): 'Maximum energy discharged'];
        end
    end
    
    %% Mathematical constructs in order to correctly measure the investment
    % Difference between investments needed.     
    for t = 1:width(xit)
        if t == 1
            Constraints = Constraints + [(yit(:,1) == xit(:,1)): 'Accumulated and New Investment equality for first year'];   % For the first year, Investment and math_indicator are the same.
        else
            Constraints = Constraints + [(yit(:,t) == xit(:,t)-xit(:,t-1)): 'Investment difference between consecutive years']; % However, for the rest, is the difference between the capacity needed in the previous year
        end                                                                                         % And this year. In case there is no new capacity added, then the incurred cost should be zero.
    end
    % Investment must not be taken back.
    % As experienced, we need to add a constraint relating to Investment that the following year, there must be at least the same amount of already invested infraestructure. The problem cannot take away investment.
    for t = 2:width(xit)
        Constraints = Constraints + [(xit(:,t) >= xit(:,t-1)): 'Investment cannot be taken back'];
    end
    % Annuities costs 
    % Candidate Thermal Generators
    for t = 1:width(GEP_Load)-1
         for i = 1:Flags.flag_thermal_cand
             rolling_window = max(1,t-Thermal_cand_data.Amortization(i)+1):t;    % This define the window of time for when the new investment is valid
             % Add actual constraint
             Constraints = Constraints + [(wit(i,t) == sum(yit(i,rolling_window))) : 'Annuities for Candidate Thermal Generators'];  % The investment is equal to the active payments on technology i
         end
    end
    % Candidate EcoThermal Generators
    for t = 1:width(GEP_Load)-1
        counter = 1;       % We keep this counter because the technology i which we are pointing in the vector Amortization is probably different than the in wit or yit
        for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
             rolling_window = max(1,t-EcoThermal_cand_data.Amortization(counter)+1):t;    % This define the window of time for when the new investment is valid
             % Add actual constraint
             Constraints = Constraints + [(wit(i,t) == sum(yit(i,rolling_window))) : 'Annuities for Candidate EcoThermal Generators'];  % The investment is equal to the active payments on technology i
             counter = counter +1;
         end
    end
    % Candidate Renewable Generators
    for t = 1:width(GEP_Load)-1
        counter = 1;       % We keep this counter because the technology i which we are pointing in the vector Amortization is probably different than the in wit or yit
        for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
             rolling_window = max(1,t-Renewables_cand_data.Amortization(counter)+1):t;    % This define the window of time for when the new investment is valid
             % Add actual constraint
             Constraints = Constraints + [(wit(i,t) == sum(yit(i,rolling_window))): 'Annuities for Candidate Renewable Generators'];  % The investment is equal to the active payments on technology i
             counter = counter +1;
         end
    end
    % Batteries
    for t = 1:width(GEP_Load)-1
        counter = 1;       % We keep this counter because the technology i which we are pointing in the vector Amortization is probably different than the in wit or yit
        for i = Flags.flag_renewables_cand+1:Flags.flag_ESS_cand
             rolling_window = max(1,t- ESS_cand_data.Amortization(counter)+1):t;    % This define the window of time for when the new investment is valid
             % Add actual constraint
             Constraints = Constraints + [(wit(i,t) == sum(yit(i,rolling_window))): 'Annuities for ESS'];  % The investment is equal to the active payments on technology i
             counter = counter +1;
         end
    end

    %% Objective Function  
    % Investment Cost
    % Information is drawn from: Capital Cost and Performance Characteristic Estimates for Utility Scale Electric Power Generating Technologies:
    % https://www.eia.gov/analysis/studies/powerplants/capitalcost/pdf/capital_cost_AEO2020.pdf
    Cand_invest_Cost = [Thermal_cand_data.InvestmentCost;EcoThermal_cand_data.InvestmentCost;Renewables_cand_data.InvestmentCost;ESS_cand_data.InvestmentCost];
    
    % Initialize Cost
    Cost = 0;
    % Add variable cost - Existing Thermal
    for t = 1:str2double(opt_year)-fst_study_horizon+1
        for i = 1:Flags.flag_thermal_exist
            Cost = Cost + Fuel_price_data.(string(Thermal_exist_data.Fuel{i}))(num2str(fst_study_horizon+t-1))*Thermal_exist_data.VarCost(i).*Exist_git(i,:,t)*((GEP_Load.LD).*8760)/((1+interest)^(fst_study_horizon+t-2021-1));
        end               % This formula reads: [AnnualVariation]*                                            *[ExistVariableCost(USD/MWh)]* *[Production(MW)]*[Duration(%)].*8760/(1+r)^t;
    end
    % Add variable cost - Existing EcoThermal           
    for t = 1:str2double(opt_year)-fst_study_horizon+1
        counter = 1;
        for i = Flags.flag_thermal_exist+1:Flags.flag_EcoThermal_exist
            if not(string(EcoThermal_exist_data.Fuel{counter}) == "Geothermal")   %% I skip the cost of Geothermal given that the renewable production cost of Geothermal is zero.
                Cost = Cost + Fuel_price_data.(string(EcoThermal_exist_data.Fuel{counter}))(num2str(fst_study_horizon+t-1))*EcoThermal_exist_data.VarCost(counter).*Exist_git(i,:,t)*((GEP_Load.LD).*8760)/((1+interest)^(fst_study_horizon+t-2021-1)); 
            end               % This formula reads: [AnnualVariation]*                                                     *[ExistVariableCost(USD/MWh)]*          *[Production(MW)]*[Duration(%)].*8760/(1+r)^t;
            counter = counter + 1;
        end
    end
    % Add variable cost - Candidate Thermal
    for t = 1:str2double(opt_year)-fst_study_horizon+1
        for i = 1:Flags.flag_thermal_cand
            Cost = Cost + Fuel_price_data.(string(Thermal_cand_data.Fuel{i}))(num2str(fst_study_horizon+t-1))*Thermal_cand_data.VarCost(i).*Cand_git(i,:,t)*((GEP_Load.LD).*8760)/((1+interest)^(fst_study_horizon+t-2021-1));
        end               % This formula reads: [AnnualVariation]*                                           *[ExistVariableCost(USD/MWh)]**[Production(MW)]*[Duration(%)].*8760/(1+r)^t;
    end
    % Add variable cost - Candidate EcoThermal
    for t = 1:str2double(opt_year)-fst_study_horizon+1
        counter = 1;
        for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
            if not(string(EcoThermal_cand_data.Fuel{counter}) == "Geothermal")        %% I skip the cost of Geothermal given that the renewable production cost of Geothermal is zero.
                Cost = Cost + Fuel_price_data.(string(EcoThermal_cand_data.Fuel{counter}))(num2str(fst_study_horizon+t-1))*EcoThermal_cand_data.VarCost(counter).*Cand_git(i,:,t)*((GEP_Load.LD).*8760)/((1+interest)^(fst_study_horizon+t-2021-1)); 
            end               % This formula reads: [AnnualVariation]*                                                    *[ExistVariableCost(USD/MWh)]*         *[Production(MW)]*[Duration(%)].*8760/(1+r)^t;
            counter = counter + 1;
        end
    end
    
    % Investment Costs
    % Add investment cost: Candidate Thermal
    for t = 1:str2double(opt_year)-fst_study_horizon+1
        for i = 1:Flags.flag_thermal_cand
            Cost = Cost + (Cand_Upper(i)*wit(i,t))*Cand_invest_Cost(i)*((interest*(1+interest)^Thermal_cand_data.Amortization(i))/((1+interest)^Thermal_cand_data.Amortization(i)-1))/((1+interest)^(fst_study_horizon+t-2021-1));
        end % Interpretation:    Cost = sum(Installed*Wit*UnitInvestmentCost*(rate(1+rate)^n)/((1+rate)^n-1)))/(1+rate)^year,       where n: is the amortization years of that machine and it might differ from t.
    end
    % Add investment cost: Candidate EcoThermal
    for t = 1:str2double(opt_year)-fst_study_horizon+1
        counter = 1;
        for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
            Cost = Cost + (Cand_Upper(i)*wit(i,t))*Cand_invest_Cost(i)*((interest*(1+interest)^EcoThermal_cand_data.Amortization(counter))/((1+interest)^EcoThermal_cand_data.Amortization(counter)-1))/((1+interest)^(fst_study_horizon+t-2021-1));
            % Interpretation:    Cost = sum(Installed*Wit*UnitInvestmentCost*(rate(1+rate)^n)/((1+rate)^n-1)))/(1+rate)^year,       where n: is the amortization years of that machine and it might differ from t.
            counter = counter + 1;
        end 
    end
    % Add investment cost: Candidate Renewables
    for t = 1:str2double(opt_year)-fst_study_horizon+1
        counter = 1;
        for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
            Cost = Cost + (Cand_Upper(i)*wit(i,t))*Cand_invest_Cost(i)*((interest*(1+interest)^Renewables_cand_data.Amortization(counter))/((1+interest)^Renewables_cand_data.Amortization(counter)-1))/((1+interest)^(fst_study_horizon+t-2021-1));
            % Interpretation:    Cost = sum(Installed*Wit*UnitInvestmentCost*(rate(1+rate)^n)/((1+rate)^n-1)))/(1+rate)^year,       where n: is the amortization years of that machine and it might differ from t.
            counter = counter + 1;
        end 
    end
    % Add investment cost: ESS
    for t = 1:str2double(opt_year)-fst_study_horizon+1
        counter = 1;
        for i = Flags.flag_renewables_cand+1:Flags.flag_ESS_cand
            Cost = Cost + (Cand_Upper(i)*wit(i,t))*Cand_invest_Cost(i)*((interest*(1+interest)^ESS_cand_data.Amortization(counter))/((1+interest)^ESS_cand_data.Amortization(counter)-1))/((1+interest)^(fst_study_horizon+t-2021-1));
            % Interpretation:    Cost = sum(Installed*Wit*UnitInvestmentCost*(rate(1+rate)^n)/((1+rate)^n-1)))/(1+rate)^year,       where n: is the amortization years of that machine and it might differ from t.
            counter = counter + 1;
        end 
    end

    % Curtailment Penalties: Both done at the same time.
    for t = 1:str2double(opt_year)-fst_study_horizon+1
        Cost = Cost + 500*sum(Curt_exist_git(:,:,t)*(8760*GEP_Load.LD))/((1+interest)^(fst_study_horizon+t-2021-1));  % FIX: The curtailment cost is fixed to 100. But, I should change it for it to be a variable. This comment applies all over the script. Especially at the UCP.
        Cost = Cost + 500*sum(Curt_cand_git(:,:,t)*(8760*GEP_Load.LD))/((1+interest)^(fst_study_horizon+t-2021-1));      
    end        % Penalty = sum(Curtailment at any given year)*PenaltyCost/(1+rate)^year
    

    %% Set options for YALMIP
    %                    verbose: amount of output information; % Our OF has big parameters so we increase Numeric Focus; % We set a time limit of 20 minutes and a Gap of 1.5%.   
    options = sdpsettings('verbose',1,'solver','gurobi','debug',1,'gurobi.NumericFocus',2,'gurobi.TimeLimit',1200,'gurobi.MIPgap',0.0015,'gurobi.Seed',rand()*1000);    % We try to add a random seed variable to check if we find alternative solutions.
    options = sdpsettings(options, 'usex0',1, 'gurobi.Heuristics',0.1);
    %                     % ,'gurobi.ImproveStartGap',0.01,'gurobi.MIPFocus',3                    
    % Let gurobi try to warm start the solution; %% Double the amount of time invested in Heuristics;% After we reach 1% optimality gap, we can look
    %                     % for better feasible solutions. % MIPFocus = 3 because the bound is not moving fast enough.
    

    % We don't solve the problem anymore, but instead create the optimizer object:
    GEP_optimizer = optimizer(Constraints,Cost,options,{sOR,ramp_capacity,Renew_Cap},{xit,Cand_git,q_Discharge_itb,q_Charge_itb,Exist_git,Cost});
end

function [Investment,Output_G_exist,Output_G_cand,Cost,GEP_optimizer] = GEP_evaluate(GEP_optimizer,sOR,ramp_capacity,DBB_info,Renew_Cap)
    import java.text.*;     v = DecimalFormat;                         % To use when printing results for easier readibility.

    % First: Evaluate the GEP Optimizer:
    [sol,errorcode,~,~,GEP_optimizer] = GEP_optimizer({sOR,ramp_capacity,Renew_Cap});
    
    % Second: call_BDD info to display better the information in Command Window and create a vector with the years
    Thermal_cand_data = DBB_info.Thermal_cand_data;
    EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
    Renewables_cand_data = DBB_info.Renewables_cand_data;
    ESS_cand_data = DBB_info.ESS_cand_data;

    Name_cand  = [Thermal_cand_data.Name;EcoThermal_cand_data.Name;Renewables_cand_data.Name;ESS_cand_data.Name];       % Candidate Names

    ColNaams = ["2021"];
    for i = 2:30
        ColNaams(i) = [convertCharsToStrings(num2str(2021-1+i))];
    end


    % Third: we will test the errorcode. In case there is no error we will convert the information to be ready to graph
    if errorcode == 0 || errorcode == 3         % Yalmip sees Gurobi running out of time as an error. We make a quickfix by including this errorcode as an admissible number.
        disp("The accumulated investment decision variable Xit: " + newline);
        disp(array2table(sol{1},'RowNames',Name_cand,'VariableNames',ColNaams));
        disp("The total cost of this GEP is: " + char(v.format(sol{6})) + " USD" + newline)   % 6: Cost 
        Cost = sol{6};
        Output_G_exist = [sol{5}];                                         % 5: Exist_git
        Output_G_cand = [sol{2};sol{3};-sol{4}];                           % 2: Cand_git; 3: q_Discharge_ih; 4: q_Charge_ih
        Investment      = sol{1};                                          % 1: xit
     
    else
        disp("The GEP case is not feasible.")
        return   % In case the mathematical formulation is not feasible, we stop the code.
    end
end

function Renew_commitment = calc_renew_commitment(opt_year)
    % Renewable quota
    % We assume that by 2020 the quota is already 50%, and by 2050 it must be
    % at least 90%, so we will do a linear function to represent this:
    % y = m*x + b
    final = 0.9; % final value commitment at the end of the study horizon
    m = (final - 0.5)/(31-1); % slope of the line
    b = 0.5 - m;              % y-intercept
    year = str2double(opt_year) - 2019;    % normalize so 2019 will be year 0
    
    Renew_commitment = m*year+ b;       % output
end

function test_days = ident_days_to_test(Global_Chronological_demand,opt_year,renewables_gen,Investment,Number,qt_days,DBB_info)  %% Warning: Might need to find a better variable name for 'Number'.
    % Find the day with the least net load: Demand - Renewable generation
    % Chronological demand of the optimized year
    Load = Global_Chronological_demand.(opt_year);

    % Retrieve renewable generation data
    Renewables_exist_data = DBB_info.Renewables_exist_data;    
    Renewables_cand_data = DBB_info.Renewables_cand_data;
    Flags = DBB_info.Flags;
    
    %% Hourly Energy Generation
    renew_gen = zeros(height(Load),1); % vector that will be used to add all of the generation of the renewables
    % Renewables generation - existing:
    for i = 1:height(Renewables_exist_data)
        renew_gen = renew_gen + renewables_gen.(Renewables_exist_data.RenewableScenario{i})*Renewables_exist_data.UpperLimit(i);
    end
    % Renewable generation - candidate:
    for i = 1:height(Renewables_cand_data)
        renew_gen = renew_gen + renewables_gen.(Renewables_cand_data.RenewableScenario{i})*Renewables_cand_data.UpperLimit(i)*Investment(Flags.flag_EcoThermal_cand+i);
    end
    
    % NetLoad
    NetLoad = Load - renew_gen;
    
    % Table for ordered information: timestamps - NetLoad
    % Use blocks flag in ordered demand base
    year = str2double(opt_year);
    % First make a first date
    t = datetime([year,1,1]);
    % Second, convert it to a number
    t = datenum(t);
    % Third, assign to time_stamps(1) as initial date
    time_stamps = NaN(Number,1);
    time_stamps(1) = t;
    % fourth, fill up the correct numbers for all the long vector time
    for ii = 2:length(time_stamps)
        time_stamps(ii) = time_stamps(ii-1) + 1/24; % for the length in time_stamps we add one hour
    end
    
    % Finally, we convert back to date format, while at the same time round to the nearest hour
    time_stamps = dateshift(datetime(time_stamps,'ConvertFrom','datenum'),'start','day','current');
    
    % Concatenate both time_stamps and NetLoad
    NetLoad_table = table(time_stamps,NetLoad);
    
    % First, we sort NetLoad_table, in ascending order, so at the begginning we will have the least NetLoad Look for the minimum 5 NetLoads
    NetLoad_table = sortrows(NetLoad_table,2);
    % Take out the first column, while at the same time apply a unique
    sort_time = unique(NetLoad_table.time_stamps,'stable');
    % Finally, take out the first qt_days unique instances:
    least_NetLoad = sort_time(1:qt_days);
    
    % Second, we sort NetLoad_table, in descending order, so at the begginng we will have the biggest NetLoad, so instances in which Load is big and renewable generation is low
    NetLoad_table = sortrows(NetLoad_table,2,'descend');
    % Take out the first column, while at the same time apply a unique, in case that we have repeating days at the begginning
    sort_time = unique(NetLoad_table.time_stamps,'stable');
    % Finally, we extract the data that we want:
    most_NetLoad = sort_time(1:qt_days);

    % Third, we will look for the instances with the most and least renewable energy production
    Renew_gen_table = table(time_stamps,renew_gen);
    % Sort the table based on renew_gen, which is the second variable:
    Renew_gen_table = sortrows(Renew_gen_table,2);
    % Take out the first column, while at the same time apply a unique:
    sort_time = unique(Renew_gen_table.time_stamps,'stable');
    % Extract the data we want:
    least_Renew_gen = sort_time(1:qt_days);
    
    % At the same time, we just reshuffle the table to order it in descending order:
    Renew_gen_table = sortrows(Renew_gen_table,2,'descend');
    % Take out the first column, while applying a unique:
    sort_time = unique(Renew_gen_table.time_stamps,'stable');
    % Extract the top qt_days data:
    most_Renew_gen = sort_time(1:qt_days);

    % Fourth, we decided we want to use Most Load. We make this change based on the logic of the four previous 
    Load_table = table(time_stamps,Load);   % So we first convert Load into a table partnered with their respective time stamps
    % Now we sort the Load table, according to Load, so according to our second variable
    Load_table = sortrows(Load_table,2,'descend');   % the descend attribute means that looking down the vector, the values descend in value. So the biggest values sit atop.
    % Now we take out the dates we desire from the first column
    sort_time = unique(Load_table.time_stamps,'stable');
    % We only extract the amount of unique days we desire according to the quantity of days we will be testing
    most_Demand = sort_time(1:qt_days);

    % Fifth, we want a set of random days.
    Random = datetime([repmat(year,qt_days,1), round((12-1).*rand(qt_days,1)) + 1,round((12-1).*rand(qt_days,1)) + 1]);


    %% Order all of the output in a single Table
    test_days = table(least_NetLoad,most_NetLoad,least_Renew_gen,most_Renew_gen,most_Demand,Random);
end

function [Week_Chrono_LNL,Week_Chrono_MNL,Week_Chrono_LRG,Week_Chrono_MRG,Week_Chrono_MDM,Week_Chrono_RAN] = Load_Chrono_to_test(opt_year,Global_Chronological_demand,test_days)
    % Extract Chronological Demand for optimized year
    Load = Global_Chronological_demand.(opt_year);
    % Table with time and load
    % First make a first date
    t = datetime([str2double(opt_year),1,1]);
    % Second, converti it to a number
    t = datenum(t);
    % Third, assign to time_stamps(1) as initial date
    time_stamps = NaN(length(Load),1);
    time_stamps(1) = t;
    % fourth, fill up the correct numbers for all the long vector time
    for ii = 2:length(time_stamps)
        time_stamps(ii) = time_stamps(ii-1) + 1/24; % for the length in time_stamps we add one hour
    end
    
    % Finally, we convert back to date format, while at the same time round
    % to the nearest hour
    time_stamps = dateshift(datetime(time_stamps,'ConvertFrom','datenum'),'start','hour','nearest');
    
    % 1. Create a matrix for the Loads of all tested days:
    % demand for Least Net Load
    Load_Chrono_LNL = zeros(24,length(test_days.least_NetLoad)); 
    % 2. Find the position for the respective days:
    for i = 1:length(test_days.least_NetLoad)
        idx = find(time_stamps == string(test_days.least_NetLoad(i)));        % the command find verifies in which position does the condition happen. In this case the condition is where does the date matches in the vector of dates
        Load_Chrono_LNL(:,i) = Load(idx:idx+23);                              % we then use that flag and extract the Load for the same position + the following 23 cases.
    end
    % 3. Interlink all of the days in a single vector to be used as Load
    Week_Chrono_LNL = [];
    for i = 1:width(Load_Chrono_LNL)
        Week_Chrono_LNL = [Week_Chrono_LNL;Load_Chrono_LNL(:,i)];
    end

    % 1. demand for Most Net Load
    Load_Chrono_MNL = zeros(24,length(test_days.most_NetLoad)); 
    % 2. Find the position for the respective days:
    for i = 1:length(test_days.most_NetLoad)
        idx = find(time_stamps == string(test_days.most_NetLoad(i)));
        Load_Chrono_MNL(:,i) = Load(idx:idx+23);
    end
     % 3. Interlink all of the days in a single vector to be used as Load
    Week_Chrono_MNL = [];
    for i = 1:width(Load_Chrono_MNL)
        Week_Chrono_MNL = [Week_Chrono_MNL;Load_Chrono_MNL(:,i)];
    end

    % 1. demand for least Renewable Generation
    Load_Chrono_LRG = zeros(24,length(test_days.least_Renew_gen)); 
    % 2. Find the position for the respective days:
    for i = 1:length(test_days.least_Renew_gen)
        idx = find(time_stamps == string(test_days.least_Renew_gen(i)));
        Load_Chrono_LRG(:,i) = Load(idx:idx+23);
    end
    % 3. Interlink all of the days in a single vector to be used as Load
    Week_Chrono_LRG = [];
    for i = 1:width(Load_Chrono_LRG)
        Week_Chrono_LRG = [Week_Chrono_LRG;Load_Chrono_LRG(:,i)];
    end

    % 1. demand for Most Renewable Generation
    Load_Chrono_MRG = zeros(24,length(test_days.most_Renew_gen)); 
    % 2. Find the position for the respective days:    
    for i = 1:length(test_days.most_Renew_gen)
        idx = find(time_stamps == string(test_days.most_Renew_gen(i)));
        Load_Chrono_MRG(:,i) = Load(idx:idx+23);
    end
    % 3. Interlink all of the days in a single vector to be used as Load
    Week_Chrono_MRG = [];
    for i = 1:width(Load_Chrono_LNL)
        Week_Chrono_MRG = [Week_Chrono_MRG;Load_Chrono_MRG(:,i)];
    end

    % 1. demand for Most Demand
    Load_Chrono_MDM = zeros(24,length(test_days.most_Demand)); 
    % 2. Find the position for the respective days:    
    for i = 1:length(test_days.most_Demand)
        idx = find(time_stamps == string(test_days.most_Demand(i)));
        Load_Chrono_MDM(:,i) = Load(idx:idx+23);
    end
    % 3. Interlink all of the days in a single vector to be used as Load
    Week_Chrono_MDM = [];
    for i = 1:width(Load_Chrono_LNL)
        Week_Chrono_MDM = [Week_Chrono_MDM;Load_Chrono_MDM(:,i)];
    end

    % 1. demand for Random
    Load_Chrono_RAN = zeros(24,length(test_days.Random)); 
    % 2. Find the position for the respective days:    
    for i = 1:length(test_days.Random)
        idx = find(time_stamps == string(test_days.Random(i)));
        Load_Chrono_RAN(:,i) = Load(idx:idx+23);
    end
    % 3. Interlink all of the days in a single vector to be used as Load
    Week_Chrono_RAN = [];
    for i = 1:width(Load_Chrono_LNL)
        Week_Chrono_RAN = [Week_Chrono_RAN;Load_Chrono_RAN(:,i)];
    end
end

function [Renewables_chrono_exist_scenario_concat,Renewables_chrono_cand_scenario_concat,Fuel_exist_index,Fuel_cand_index] = UCP_renew_gen_and_fuel(renewables_gen,dates,opt_year,DBB_info)
    % Call information to know the 
    %% Call Input Information about existing and candidate generators
    Thermal_exist_data = DBB_info.Thermal_exist_data;
    EcoThermal_exist_data = DBB_info.EcoThermal_exist_data;
    Renewables_exist_data = DBB_info.Renewables_exist_data;    
    Thermal_cand_data = DBB_info.Thermal_cand_data;
    EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
    Renewables_cand_data = DBB_info.Renewables_cand_data;
    Fuel_price_data = DBB_info.Fuel_price_data;
    Flags = DBB_info.Flags;

    %%% RENEWABLES
    %% 1. Extract the information for the generators
    % First make a first date
    t = datetime([str2double(opt_year),1,1]);
    % Second, converti it to a number
    t = datenum(t);
    % Third, assign to time_stamps(1) as initial date
    time_stamps = NaN(height(renewables_gen),1);
    time_stamps(1) = t;
    % fourth, fill up the correct numbers for all the long vector time
    for ii = 2:length(time_stamps)
        time_stamps(ii) = time_stamps(ii-1) + 1/24; % for the length in time_stamps we add one hour
    end
    % Finally, we convert back to date format, while at the same time round
    % to the nearest hour
    time_stamps = dateshift(datetime(time_stamps,'ConvertFrom','datenum'),'start','hour','nearest');

    % We now delete the variable time stamps from renewables_gen, making a
    % new variable with a slightly different name
    Renewables_gen  = removevars(renewables_gen,{'time_stamp'});
    % Immediately after, insert the new 'time_stamp' column:
    Renewables_gen  = addvars(Renewables_gen,time_stamps,'Before','block_flag');
    
    % Extract the ids of the dates that we need to extract
    id   = zeros(length(dates),1);
    for i = 1:length(dates)                                     % Warning: no need to fix when stumbling into trouble.          
        id(i) = find(Renewables_gen.time_stamps == dates(i));   % Warning: it stumbles into trouble when used outside the main code because the test daays change between loop,
    end                                                         % but when tested, the variable test_days is not updated while Renewables_gen is.

    %% 2. Extract the Chronological Generation for the Existing generators
    % Solar Chronological resource
    Renewables_chrono_exist_scenario = nan(24,height(Renewables_exist_data),length(dates));  % Make a NAN 3-D array to be replaced with information extracted from the chronological renewable scenario
    for j = 1:length(id) 
        counter = 1;        % Given that we need to start the counter at the top line of the RenewablesExist.xlsb file
        for i = Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist           % for depending on the amount of existing solar generators
            Renewables_exist = Renewables_gen.(Renewables_exist_data.RenewableScenario{counter});     % We extract the information for the specific renewable scenario (all year)
            Renewables_exist = Renewables_exist(id(j):id(j)+23);                                      % We cut the information just for the dates we are interested in
            Renewables_chrono_exist_scenario(:,counter,j) = Renewables_exist;                         % We add to the master array
            counter = counter + 1;
        end
    end
    % Make the 3-D array into a 2-D array by concatenating the matrices into a single matrix.
    Renewables_chrono_exist_scenario_concat = []; % Empty matrix variable that in which the other information will be pasted to.
    for i = 1:size(Renewables_chrono_exist_scenario,3)
        Renewables_chrono_exist_scenario_concat = [Renewables_chrono_exist_scenario_concat;Renewables_chrono_exist_scenario(:,:,i)];  %% We now just append below the information of all existing scenarios in a single 2-D Matrix
    end

    %% 3. Extract the Chronological Generation for the Candidate generators
    % Solar Chronological resource
    Renewables_chrono_cand_scenario = nan(24,height(Renewables_cand_data),length(dates));  % Make a NAN 3-D array to be replaced with information extracted from the chronological renewable scenario
    for j = 1:length(id) 
        counter = 1;        % Given that we need to start the counter at the top line of the SolarExit.xlsb file
        for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand          % for depending on the amount of existing solar generators
            Renewables_cand = Renewables_gen.(Renewables_cand_data.RenewableScenario{counter});     % We extract the information for the specific renewable scenario (all year)
            Renewables_cand = Renewables_cand(id(j):id(j)+23);                                      % We cut the information just for the dates we are interested in
            Renewables_chrono_cand_scenario(:,counter,j) = Renewables_cand;                         % We add to the master array
            counter = counter + 1;
        end
    end
    % Make the 3-D array into a 2-D array by concatenating the matrices into a single matrix.
    Renewables_chrono_cand_scenario_concat = []; % Empty matrix variable that in which the other information will be pasted to.
    for i = 1:size(Renewables_chrono_exist_scenario,3)
        Renewables_chrono_cand_scenario_concat = [Renewables_chrono_cand_scenario_concat;Renewables_chrono_cand_scenario(:,:,i)];  %% We now just append below the information of all existing scenarios in a single 2-D Matrix
    end

    %%% FUEL PRICES
    % Existing Fuel Forecast: Thermal
    Fuel_exist_index = [];
    for i = 1:Flags.flag_thermal_exist             % First: finds the fuel of technology i. Then, extracts the Fuel_price_data for that fuel.
        Fuel_exist_index = [Fuel_exist_index,Fuel_price_data.(string(Thermal_exist_data.Fuel{i}))(opt_year)]; % Finally, it obtains just a number indexed for year opt_year
    end
    % Existing Fuel Variation: EcoThermal
    for i = 1: Flags.flag_EcoThermal_exist - Flags.flag_thermal_exist
        Fuel_exist_index = [Fuel_exist_index,Fuel_price_data.(string(EcoThermal_exist_data.Fuel{i}))(opt_year)];
    end
    
    % Candidate Fuel Variation: Thermal
    Fuel_cand_index = [];
    for i = 1:Flags.flag_thermal_cand
        Fuel_cand_index = [Fuel_cand_index,Fuel_price_data.(string(Thermal_cand_data.Fuel{i}))(opt_year)];
    end
    % Candidate Fuel Variatin: EcoThermal
    for i = 1: Flags.flag_EcoThermal_cand - Flags.flag_thermal_cand
        Fuel_cand_index = [Fuel_cand_index,Fuel_price_data.(string(EcoThermal_cand_data.Fuel{i}))(opt_year)];
    end
end

function UCP_optimizer = UCP_optimizer_creator(size_Load,DBB_info)
    % Nomenclature:
    % i: Generators. Will be used both for aggregated generators in Existing or Candidate decision variable object.
    % h: hours.
    % t: years. -> However, there is not instance in here where we iterate over years. We only use the optimized year variable to keep track of the year in which we are right now.
    
    % Call BDD information/data
    Thermal_exist_data = DBB_info.Thermal_exist_data;
    EcoThermal_exist_data = DBB_info.EcoThermal_exist_data;
    Renewables_exist_data = DBB_info.Renewables_exist_data;    
    Thermal_cand_data = DBB_info.Thermal_cand_data;
    EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
    Renewables_cand_data = DBB_info.Renewables_cand_data;
    ESS_cand_data = DBB_info.ESS_cand_data;
    Flags = DBB_info.Flags;

    % General Information
    % Existing Generator Bounds
    Exist_Upper = [Thermal_exist_data.UpperLimit; EcoThermal_exist_data.UpperLimit; Renewables_exist_data.UpperLimit]; % Upper Generation Bounds Existing
    Exist_Lower = [Thermal_exist_data.LowerLimit; EcoThermal_exist_data.LowerLimit; zeros(height(Renewables_exist_data),1)]; % Lower Generation Bounds Existing
    
    % Candidate Generator Bounds
    Cand_Upper = [Thermal_cand_data.UpperLimit; EcoThermal_cand_data.UpperLimit; Renewables_cand_data.UpperLimit]; % Upper Generation Bounds Candidate
    Cand_Lower = [Thermal_cand_data.LowerLimit; EcoThermal_cand_data.LowerLimit; zeros(height(Renewables_cand_data),1)]; % Lower Generation Bounds Candidate
      
    %% Mathematical Formulation
    % Decision variables
    Exist_git       = sdpvar(Flags.quant_exist,size_Load,'full');            % continuous variable to define generation output from exist generation
    Cand_git        = sdpvar(Flags.flag_renewables_cand,size_Load,'full');   % continuous variable to define generation output from future generation, excluding ESS
    CommitExist_zit = binvar(Flags.quant_exist,size_Load,'full');            % binary variable to determine if a Existing unit is committed or not at a given moment
    CommitCand_zit  = intvar(Flags.flag_renewables_cand,size_Load,'full');   % integer variable to determine if a Candidate unit is committed or not at a given moment, excluding ESS
    StartExist_yit  = binvar(Flags.quant_exist,size_Load,'full');            % binary variable to flag if an Existing Thermal Technology starts-up
    StartCand_yit   = intvar(Flags.flag_renewables_cand,size_Load,'full');   % integer variable to flag if a Candidate Thermal Technology starts-up, excluding ESS
    ShutD_exist_xit = binvar(Flags.quant_exist,size_Load,'full');            % binary variable to flag if an Existing Thermal Technology Shuts-down
    ShutD_cand_xit  = intvar(Flags.flag_renewables_cand,size_Load,'full');   % integer variable to flag if a Candidate Thermal Technology Shuts-down, excluding ESS
    q_SOC_ih        = sdpvar(height(ESS_cand_data),size_Load,'full');        % continuous variable decision of amount of energy stored in each type of battery
    q_Discharge_ih  = sdpvar(height(ESS_cand_data),size_Load,'full');        % continuous variable of Discharge by ESS i at hour h
    q_Charge_ih     = sdpvar(height(ESS_cand_data),size_Load,'full');        % continuous variable of Charge by ESS i at hour h
    % Curtailment decision variables
    Curt_Exist_git  = sdpvar(Flags.flag_renewables_exist-Flags.flag_EcoThermal_exist,size_Load,'full');    % continuous variable to keep track of curtailment of existing renewables
    Curt_Cand_git  = sdpvar(Flags.flag_renewables_cand-Flags.flag_EcoThermal_cand,size_Load,'full');       % continuous variable to keep track of curtailment of existing renewables

    % Parameters
    % Generators
    Investment = sdpvar(Flags.quant_cand,1,'full'); % Parameter for the vector of Investment passed by the GEP
    % Load
    Load       = sdpvar(size_Load,1,'full');        % Parameter for load changing per UCP tested. This changes by year, and by iteration of the GEP
    % Renewables
    Renewables_chrono_exist_scenario = sdpvar(size_Load,Flags.flag_renewables_exist-Flags.flag_EcoThermal_exist,'full');    % Parameter that will dictated the resource available for that specific hour for existing renewables
    Renewables_chrono_cand_scenario  = sdpvar(size_Load,Flags.flag_renewables_cand-Flags.flag_EcoThermal_cand,'full');      % Parameter that will dictated the resource available for that specific hour for candidate renewables
    % Fuel Price
    Fuel_exist_index = sdpvar(1,Flags.flag_EcoThermal_exist,'full');     % Parameter that will be multiplied against the energy generated by existing thermal power, including ecothermal.
    Fuel_cand_index  = sdpvar(1,Flags.flag_EcoThermal_cand,'full');      % Parameter that will be multiplied against the energy generated by candidate thermal power, including ecothermal.
    
    %% Constraints
    % General Constraints
    Constraints = [];
    % Power Balance Constraint
    for h = 1:size_Load
        Constraints = Constraints + [(sum(Exist_git(:,h)) + sum(Cand_git(:,h)) + sum(q_Discharge_ih(:,h)) == Load(h) + sum(q_Charge_ih(:,h))) : 'Power Balance Constraint']; % 
    end

    % Generator Bounds: Min and Max Constraints
    % Existing Generator Bounds
    for h = 1:size_Load          % This constraint reads: [Binvar*LowerBound <= Generation <= Binvar*UpperBound ]; In other words, if a technology is committed its generation must be within bounds. If not, it must be zero.
        Constraints = Constraints + [(CommitExist_zit(:,h).*Exist_Lower <= Exist_git(:,h) <= CommitExist_zit(:,h).*Exist_Upper): 'Existing Generator Bounds'];
    end
    % Candidate Generator Bounds
    for h = 1:size_Load          % This constraint reads: [Binvar*LowerBound <= Generation <= Binvar*UpperBound ]; In other words, if a technology is committed its generation must be within bounds. If not, it must be zero.
        Constraints = Constraints + [(CommitCand_zit(:,h).*Cand_Lower <= Cand_git(:,h) <= CommitCand_zit(:,h).*Cand_Upper): 'Candidate Generator Bounds']; 
    end

    % Commitment Integer Upper bound
    for h = 1:size_Load          % This constraints reads: zit <= Ni, wher zit is the commitment and Ni is the installed number of technology i
        Constraints = Constraints + [(CommitCand_zit(:,h) <= Investment(1:Flags.flag_renewables_cand)) : 'Commitment at most to Installed Candidate'];
    end


    %% Thermal Constraints
    % Ramp up and Ramp down Constraints
    % Existing Thermal Generators
    for i = 1:Flags.flag_thermal_exist
        for h = 2:size_Load
            Constraints = Constraints + [(Exist_git(i,h-1) - Thermal_exist_data.RampDown(i)*CommitExist_zit(i,h) <= Exist_git(i,h) <= Exist_git(i,h-1) + Thermal_exist_data.RampUp(i)*CommitExist_zit(i,h)): 'Existing Thermal Ramp up and Down Constraint'];   % Ramp-up and ramp-down constraint
        end
    end
    % Existing EcoThermal Generators
    counter = 1;        % counter to distinguish when counting for the Ecothermal data and the Git matrix.
    for i = Flags.flag_thermal_exist+1:Flags.flag_EcoThermal_exist
        for h = 2:size_Load
            Constraints = Constraints + [(Exist_git(i,h-1) - EcoThermal_exist_data.RampDown(counter)*CommitExist_zit(i,h) <= Exist_git(i,h) <= Exist_git(i,h-1) + EcoThermal_exist_data.RampUp(counter)*CommitExist_zit(i,h)): 'Existing EcoThermal Ramp up and down Constraint'];   % Ramp-up and ramp-down constraint
        end
        counter = counter + 1;
    end
    % Cand Thermal Generators
    for i = 1:Flags.flag_thermal_cand
        for h = 2:size_Load          % Ramp-up and ramp-down constraint
            Constraints = Constraints + [(Cand_git(i,h-1) - Thermal_cand_data.RampDown(i)*CommitCand_zit(i,h) <= Cand_git(i,h) <= Cand_git(i,h-1) + Thermal_cand_data.RampUp(i)*CommitCand_zit(i,h)) : 'Candidate Thermal Ramp up and down Constraint'];
        end
    end
    % Candidate EcoThermal Generators
    counter = 1;        % counter to distinguish when counting for the Ecothermal data and the Git matrix.
    for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
        for h = 2:size_Load          % Ramp-up and ramp-down constraint
            Constraints = Constraints + [(Cand_git(i,h-1) - EcoThermal_cand_data.RampDown(counter)*CommitCand_zit(i,h) <= Cand_git(i,h) <= Cand_git(i,h-1) + EcoThermal_cand_data.RampUp(counter)*CommitCand_zit(i,h)) :'Candidate EcoThermal Ramp up and down Constraint'];
        end
        counter = counter + 1;
    end

    % Start Indicator: All Technologies at the same time.          % Observation: We are defining a flag for start for Renewables, but we will defined a Cost Zero, so the solver might take whatever value it decides to be okay.
    for h = 2:size_Load
        Constraints = Constraints + [(CommitExist_zit(:,h) <= CommitExist_zit(:,h-1) + StartExist_yit(:,h)) : 'Startup Indicator Existing Generators'];
        Constraints = Constraints + [(CommitCand_zit(:,h) <=  CommitCand_zit(:,h-1)  + StartCand_yit(:,h)): 'Startup Indicator Candidate Generators'];
    end
    % Shutdown Indicator: All Techonologies done at the same time. % Observation: We are defining a flag for start for Renewables, but we will defined a Cost Zero, so the solver might take whatever value it decides to be okay.                             
    for h = 2:size_Load                                
        Constraints = Constraints + [(CommitExist_zit(:,h) >= CommitExist_zit(:,h-1) - ShutD_exist_xit(:,h)): 'Shutdown Indicator Existing Generators'];
        Constraints = Constraints + [(CommitCand_zit(:,h) >= CommitCand_zit(:,h-1) - ShutD_cand_xit(:,h)): 'Shutdown Indicator Candidate Generators'];
    end    

    % MinUp and MinDown time
    % MinUp and MinDown Existing Thermal Generators
    % Existing Thermal Generators
    for h = 2:size_Load
        for i = 1:Flags.flag_thermal_exist
            rangeUp   = h:min(size_Load,h+Thermal_exist_data.MinUp(i)-1);
            rangeDown = h:min(size_Load,h+Thermal_exist_data.MinDown(i)-1);
            % Add actual Constraints
            Constraints = Constraints + [(CommitExist_zit(i,rangeUp) >= StartExist_yit(i,h)): 'MinUp Time Existing Thermal'];        % MinUp time
            Constraints = Constraints + [(CommitExist_zit(i,rangeDown) <= 1 - ShutD_exist_xit(i,h)): 'MinDown Time Existing Thermal']; % MinDown time
        end
    end
    % Existing EcoThermal Generators
    for h = 2:size_Load
        counter = 1;        % counter to distinguish EcoThermal and Git.
        for i = Flags.flag_thermal_exist+1:Flags.flag_EcoThermal_exist
            rangeUp =   h:min(size_Load,h+EcoThermal_exist_data.MinUp(counter)-1);
            rangeDown = h:min(size_Load,h+EcoThermal_exist_data.MinDown(counter)-1);
            % Add actual Constraints
            Constraints = Constraints + [(CommitExist_zit(i,rangeUp) >= StartExist_yit(i,h)): 'MinUp Time Existing EcoThermal'];        % MinUp time
            Constraints = Constraints + [(CommitExist_zit(i,rangeDown) <= 1 - ShutD_exist_xit(i,h)): 'MinDown Time Existing EcoThermal']; % MinDown time
            counter = counter + 1;
        end
    end
    % Candidate Thermal Generators
    for h = 2:size_Load
        for i = 1:Flags.flag_thermal_cand
            rangeUp =   h:min(size_Load,h+Thermal_cand_data.MinUp(i)-1);
            rangeDown = h:min(size_Load,h+Thermal_cand_data.MinDown(i)-1);
            % Add actual Constraints
            Constraints = Constraints + [(CommitCand_zit(i,rangeUp) >= StartCand_yit(i,h)): 'MinUp Time Candidate Thermal'];         % MinUp time
            Constraints = Constraints + [(CommitCand_zit(i,rangeDown) <= Investment(i) - ShutD_cand_xit(i,h)): 'MinDown Time Candidate EcoThermal'];  % MinDown time
        end
    end
    % Candidate EcoThermal Generators
    for h = 2:size_Load
        counter = 1;        % counter to distinguish EcoThermal and Git.
        for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
            rangeUp =   h:min(size_Load,h+EcoThermal_cand_data.MinUp(counter)-1);
            rangeDown = h:min(size_Load,h+EcoThermal_cand_data.MinDown(counter)-1);
            % Add actual Constraints
            Constraints = Constraints + [(CommitCand_zit(i,rangeUp) >= StartCand_yit(i,h)): 'MinUp Time Candidate EcoThermal'];        % MinUp time
            Constraints = Constraints + [(CommitCand_zit(i,rangeDown) <= Investment(i) - ShutD_cand_xit(i,h)): 'MinDown Time Candidate EcoThermal']; % MinDown time
            counter = counter + 1;
        end
    end
    
    %% Renewables generation resource: This constraint is read that the generation from renewables needs to be at most the available resource
    % Existing renewables
    counter = 1; % fictitious variable to start the count for the length of renewables
    for i = Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist
        Constraints = Constraints + [(Exist_git(i,:)' <= Exist_Upper(i).*Renewables_chrono_exist_scenario(:,counter)): 'Existing Renewable Available Resource'];
        counter = counter +1;
    end
    % Candidate renewables
    counter = 1;
    for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
        Constraints = Constraints + [(Cand_git(i,:)' <= (Investment(i)*Cand_Upper(i)).*Renewables_chrono_cand_scenario(:,counter)): 'Candidate Renewable Available Resource'];
        counter = counter +1;
    end

    % Constraints to keep track of the curtailment
    % Existing renewables
    counter = 1;
    for i = Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist
        Constraints = Constraints + [(Curt_Exist_git(counter,:)' == Exist_Upper(i).*Renewables_chrono_exist_scenario(:,counter) - Exist_git(i,:)'): 'Existing Renewable Curtailment'];
        counter = counter +1;       % Given that Curt_Exist_git dimensions are only for renewables, we keep track with 'counter' instead of 'i'
    end
    % Candidate renewables
    counter = 1;
    for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
        Constraints = Constraints + [(Curt_Cand_git(counter,:)' == (Investment(i)*Cand_Upper(i)).*Renewables_chrono_cand_scenario(:,counter) - Cand_git(i,:)'): 'Candidate Renewable Curtailment'];
        counter = counter +1;       % Given that Curt_Exist_git dimensions are only for renewables, we keep track with 'counter' instead of 'i'
    end                

    %% ESS Additional constraints
    % State of charge constraints
    for i = 1:height(ESS_cand_data)
        for h = 2:size_Load
            Constraints = Constraints + [(q_SOC_ih(i,h) == q_SOC_ih(i,h-1) - q_Discharge_ih(i,h) + ESS_cand_data.Efficiency(i)*q_Charge_ih(i,h)): 'SOC energy conservation'];
        end
    end
    % Initial and final state of charge per day
    hours = [0:24:size_Load];
    hours(1)= 1;
    counter = Flags.flag_renewables_cand+1;
    for i = 1:height(ESS_cand_data)
        for h = 1:length(hours)
            Constraints = Constraints + [(q_SOC_ih(i,hours(h)) == 1/2*ESS_cand_data.Duration(i)*ESS_cand_data.UpperLimit(i)*Investment(counter)): 'SOC same at the end of everyday'];
        end
        counter = counter +1;
    end
    % Maximum energy stored at the batteries
    counter = Flags.flag_renewables_cand+1;
    for i = 1:height(ESS_cand_data)
        for h = 1:size_Load
            Constraints = Constraints + [(0 <= q_SOC_ih(i,h) <= ESS_cand_data.Duration(i)*ESS_cand_data.UpperLimit(i)*Investment(counter)): 'Energy Bounds in Batteries'];
        end
        counter = counter +1;
    end
    % Charge and Discharge bounds
    counter = Flags.flag_renewables_cand+1;
    for i = 1:height(ESS_cand_data)
        for h = 1:size_Load
            Constraints = Constraints + [(0 <= q_Discharge_ih(i,h) <= ESS_cand_data.UpperLimit(i)*Investment(counter)): 'Batteries Discharge Power Bound'];
            Constraints = Constraints + [(0 <= q_Charge_ih(i,h) <= ESS_cand_data.UpperLimit(i)*Investment(counter)): 'Batteries Charge Power Bound' ];
        end
        counter = counter +1;
    end

    %% Security Constraints
    % Maximum Curtailment: um(Curt_Exist_git) + sum(Curt_Cand_exist) <= 0.6(Exist_ResourceAvailability*RenewablesInstalledCapacity + Cand_RenewablesInstalledCapacity*Investment*ResourceAvailability)
    Constraints = Constraints + [sum(sum(Curt_Exist_git)) + sum(sum(Curt_Cand_git)) <= 0.4*(sum(Renewables_chrono_exist_scenario*Exist_Upper(Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist)) + sum(Renewables_chrono_cand_scenario*(Investment(Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand).*Cand_Upper(Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand))))];

    %% Cap for UCP Payment
    % We will define a maximum payment for UCP within a week. The cap is 1billion USD
    % This Constraint was built after the Objective Function. However, it will be based on the structure of the Objective Function.
    % We define the LHS as the sum of costs, while the RHS is just the cap = 1billion USD
    Cap_LHS = 0;
    % Add Thermal generation cost - Existing
    for i = 1:Flags.flag_thermal_exist
        Cap_LHS = Cap_LHS + sum(Fuel_exist_index(i)*Thermal_exist_data.VarCost(i)*Exist_git(i,:)); % This formula reads: sum([VariableCost(USD/MWh)]*[Production(MW)]);
    end
    % Add EcoThermal generation cost - Existing
    counter = 1;
    for i = Flags.flag_thermal_exist+1:Flags.flag_EcoThermal_exist
        Cap_LHS = Cap_LHS + sum(Fuel_exist_index(i)*EcoThermal_exist_data.VarCost(counter)*Exist_git(i,:)); % This formula reads: sum([VariableCost(USD/MWh)]*[Production(MW)]);
        counter = counter + 1;
    end
    % Add Thermal generation cost - Candidate
    for i = 1:Flags.flag_thermal_cand
        Cap_LHS = Cap_LHS + sum(Fuel_cand_index(i)*Thermal_cand_data.VarCost(i)*Cand_git(i,:)); % This formula reads: sum([VariableCost(USD/MWh)]*[Production(MW)]);
    end
    % Add EcoThermal generation cost - Candidate
    counter = 1;
    for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
        Cap_LHS = Cap_LHS + sum(Fuel_cand_index(i)*EcoThermal_cand_data.VarCost(counter)*Cand_git(i,:)); % This formula reads: sum([VariableCost(USD/MWh)]*[Production(MW)]);
        counter = counter + 1;
    end
    % On the section "Add generation cost" the second term is a vector, so we will insert a sum function to add the whole vector and minimize a single Objective Function
    Cap_LHS = sum(Cap_LHS);
    
    % Add StartUp Cost and Shutdown Cost                % We only sum from the second hour onwards, because for the first hour it might not be defined.
    % Existing Thermal Generators
    for i = 1:Flags.flag_thermal_exist
        Cap_LHS = Cap_LHS + sum(StartExist_yit(i,2:end))*Thermal_exist_data.StartUpCost(i);
        Cap_LHS = Cap_LHS + sum(ShutD_exist_xit(i,2:end))*Thermal_exist_data.StartUpCost(i);
    end
    % Existing EcoThermal Generators
    counter = 1;
    for i = Flags.flag_thermal_exist+1:Flags.flag_EcoThermal_exist
        Cap_LHS = Cap_LHS + sum(StartExist_yit(i,2:end))*EcoThermal_exist_data.StartUpCost(counter);
        Cap_LHS = Cap_LHS + sum(ShutD_exist_xit(i,2:end))*EcoThermal_exist_data.StartUpCost(counter);
        counter = counter +1;
    end
    % Candidate Thermal Generators
    for i = 1:Flags.flag_thermal_cand
        Cap_LHS = Cap_LHS + sum(StartCand_yit(i,2:end))*Thermal_cand_data.StartUpCost(i);
        Cap_LHS = Cap_LHS + sum(ShutD_cand_xit(i,2:end))*Thermal_cand_data.StartUpCost(i);
    end
    % Candidate EcoThermal Generators
    counter = 1;
    for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
        Cap_LHS = Cap_LHS + sum(StartCand_yit(i,2:end))*EcoThermal_cand_data.StartUpCost(counter);
        Cap_LHS = Cap_LHS + sum(ShutD_cand_xit(i,2:end))*EcoThermal_cand_data.StartUpCost(counter);
        counter = counter +1;
    end
    % Add Curtailment penalty cost: 1. When we first sum(Curtailment) we obtain a vector. 2. When we sum(sum(Curtailment)) we get a number. This number is added to the cost multiplied by 100 USD/MWh.
    % Existing Curtailment
    Cap_LHS = Cap_LHS + 300*sum(sum(Curt_Exist_git));
    % Candidate Curtailment
    Cap_LHS = Cap_LHS + 300*sum(sum(Curt_Cand_git));

    % Finally add the constraint
    Constraints = Constraints + [(Cap_LHS <= 1000000000): 'Cap UCP Financial Resource'];

    %% Objective Function    % FIX: THE UCP IS NOT INDEXING THE COST OF FUEL ACROSS TIME. Refer to v8. UCP_optimizer(sOR,ramp,opt_year). This also affects the 'Cap UCP Financial Resource' contraint.
    % Initialize Cost
    Cost = 0;
    % Add Thermal generation cost - Existing
    for i = 1:Flags.flag_thermal_exist
        Cost = Cost + sum(Fuel_exist_index(i)*Thermal_exist_data.VarCost(i)*Exist_git(i,:)); % This formula reads: sum([VariableCost(USD/MWh)]*[Production(MW)]);
    end
    % Add EcoThermal generation cost - Existing
    counter = 1;
    for i = Flags.flag_thermal_exist+1:Flags.flag_EcoThermal_exist
        Cost = Cost + sum(Fuel_exist_index(i)*EcoThermal_exist_data.VarCost(counter)*Exist_git(i,:)); % This formula reads: sum([VariableCost(USD/MWh)]*[Production(MW)]);
        counter = counter + 1;
    end
    % Add Thermal generation cost - Candidate
    for i = 1:Flags.flag_thermal_cand
        Cost = Cost + sum(Fuel_cand_index(i)*Thermal_cand_data.VarCost(i)*Cand_git(i,:)); % This formula reads: sum([VariableCost(USD/MWh)]*[Production(MW)]);
    end
    % Add EcoThermal generation cost - Candidate
    counter = 1;
    for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
        Cost = Cost + sum(Fuel_cand_index(i)*EcoThermal_cand_data.VarCost(counter)*Cand_git(i,:)); % This formula reads: sum([VariableCost(USD/MWh)]*[Production(MW)]);
        counter = counter + 1;
    end
    % On the section "Add generation cost" the second term is a vector, so we will insert a sum function to add the whole vector and minimize a single Objective Function
    Cost = sum(Cost);
    
    % Add StartUp Cost and Shutdown Cost                % We only sum from the second hour onwards, because for the first hour it might not be defined.
    % Existing Thermal Generators
    for i = 1:Flags.flag_thermal_exist
        Cost = Cost + sum(StartExist_yit(i,2:end))*Thermal_exist_data.StartUpCost(i);
        Cost = Cost + sum(ShutD_exist_xit(i,2:end))*Thermal_exist_data.StartUpCost(i);
    end
    % Existing EcoThermal Generators
    counter = 1;
    for i = Flags.flag_thermal_exist+1:Flags.flag_EcoThermal_exist
        Cost = Cost + sum(StartExist_yit(i,2:end))*EcoThermal_exist_data.StartUpCost(counter);
        Cost = Cost + sum(ShutD_exist_xit(i,2:end))*EcoThermal_exist_data.StartUpCost(counter);
        counter = counter +1;
    end
    % Candidate Thermal Generators
    for i = 1:Flags.flag_thermal_cand
        Cost = Cost + sum(StartCand_yit(i,2:end))*Thermal_cand_data.StartUpCost(i);
        Cost = Cost + sum(ShutD_cand_xit(i,2:end))*Thermal_cand_data.StartUpCost(i);
    end
    % Candidate EcoThermal Generators
    counter = 1;
    for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
        Cost = Cost + sum(StartCand_yit(i,2:end))*EcoThermal_cand_data.StartUpCost(counter);
        Cost = Cost + sum(ShutD_cand_xit(i,2:end))*EcoThermal_cand_data.StartUpCost(counter);
        counter = counter +1;
    end
    % Add Curtailment penalty cost: 1. When we first sum(Curtailment) we obtain a vector. 2. When we sum(sum(Curtailment)) we get a number. This number is added to the cost multiplied by 100 USD/MWh.
    % Existing Curtailment
    Cost = Cost + 300*sum(sum(Curt_Exist_git));
    % Candidate Curtailment
    Cost = Cost + 300*sum(sum(Curt_Cand_git));

    %% Solve the problem
    % Set options for YALMIP
    options = sdpsettings('verbose',1,'solver','gurobi','debug',1,'gurobi.NumericFocus',1,'gurobi.MIPGap',0.003,'gurobi.TimeLimit',1200,'gurobi.Seed',rand()*1000,'gurobi.Heuristics',0.1,'usex0',1);

    % Create optimizer object
    UCP_optimizer = optimizer(Constraints,[],options,{Investment,Load,Renewables_chrono_exist_scenario,Renewables_chrono_cand_scenario,Fuel_exist_index,Fuel_cand_index},{Exist_git,Cand_git,q_Discharge_ih,q_Charge_ih,Cost,Curt_Exist_git,Curt_Cand_git});
end

function [Output_G_cand,Output_G_exist,errorcode,UCP_optimizer,Case_Curtailment] = UCP_evaluate(renewables_gen,dates,opt_year,Investment,Load,UCP_optimizer,UCP_type,DBB_info)
    import java.text.*;     v = DecimalFormat;                         % To use when printing results for easier readibility.
    
    % First: we calculate the renewables and fuel variation to use with the UCP optimizer
    % For that we will call another function that does that
    [Renewables_chrono_exist_scenario,Renewables_chrono_cand_scenario,Fuel_Forecast_exist,Fuel_Forecast_cand] = UCP_renew_gen_and_fuel(renewables_gen,dates,opt_year,DBB_info);
    
    % Second: With all of this information and the Load and investment that are inputs from the original function, we can use the optimizer.
    [sol,errorcode,~,~,UCP_optimizer] = UCP_optimizer({Investment,Load,Renewables_chrono_exist_scenario,Renewables_chrono_cand_scenario,Fuel_Forecast_exist,Fuel_Forecast_cand});

    % Third: we will test the errorcode. In case there is no error we will convert the information to be ready to graph
    if errorcode == 0 || errorcode == 3
        disp("The cost associated with this UCP (" + UCP_type +") is: " + char(v.format(sol{5})) + " USD")   % 5: Cost 
        Output_G_exist = [sol{1}];                                         % 1: Exist_git
        Output_G_cand = [sol{2};sol{3};-sol{4}];                           % 2: Cand_git; 3: q_Discharge_ih; 4: q_Charge_ih
        % Curtailment output to reach a decision about how to give back information to the GEP:
        Thermal_exist_data = DBB_info.Thermal_exist_data;
        EcoThermal_exist_data = DBB_info.EcoThermal_exist_data;
        Renewables_exist_data = DBB_info.Renewables_exist_data;    
        Thermal_cand_data = DBB_info.Thermal_cand_data;
        EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
        Renewables_cand_data = DBB_info.Renewables_cand_data;
        Flags = DBB_info.Flags;
    
        % General Information
        % Existing Generator Bounds
        Exist_Upper = [Thermal_exist_data.UpperLimit; EcoThermal_exist_data.UpperLimit; Renewables_exist_data.UpperLimit]; % Upper Generation Bounds Existing
        
        % Candidate Generator Bounds
        Cand_Upper = [Thermal_cand_data.UpperLimit; EcoThermal_cand_data.UpperLimit; Renewables_cand_data.UpperLimit]; % Upper Generation Bounds Candidate
        Case_Curtailment = (sum(sum(sol{6})) + sum(sum(sol{7})))/(sum(Renewables_chrono_exist_scenario*Exist_Upper(Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist)) + sum(Renewables_chrono_cand_scenario*(Investment(Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand).*Cand_Upper(Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand))));
        if Case_Curtailment >= 0.2
            disp("This case curtailment is: " + string(v.format((sum(sum(sol{6})) + sum(sum(sol{7})))/(sum(Renewables_chrono_exist_scenario*Exist_Upper(Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist)) + sum(Renewables_chrono_cand_scenario*(Investment(Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand).*Cand_Upper(Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand)))))))
        end
        
    else
        disp("This case is not feasible.")
        Output_G_cand = nan;
        Output_G_exist = nan;
        Case_Curtailment = nan;
    end
end

function fig = graph_EnergyProduction(Output_G_exist,Output_G_cand,opt_year,graph_name_flag,sOR,eval_criteria,visibility,GEP_Load,DBB_info,ramp_capacity)
    import java.text.*;     v = DecimalFormat;  % To use when printing results for easier readibility.
    %% Extract information of generators Name
    % Call Input Information about existing and candidate generators
    Thermal_exist_data = DBB_info.Thermal_exist_data;
    EcoThermal_exist_data = DBB_info.EcoThermal_exist_data;
    Renewables_exist_data = DBB_info.Renewables_exist_data;    
    Thermal_cand_data = DBB_info.Thermal_cand_data;
    EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
    Renewables_cand_data = DBB_info.Renewables_cand_data;
    ESS_cand_data = DBB_info.ESS_cand_data;
    Flags = DBB_info.Flags;

    % Name of Generators: Existing
    Name_exist = [Thermal_exist_data.Name;EcoThermal_exist_data.Name;Renewables_exist_data.Name];
    % Name of Generators: Candidate
    Name_cand  = [Thermal_cand_data.Name;EcoThermal_cand_data.Name;Renewables_cand_data.Name;ESS_cand_data.Name + " discharges"];       % We add to the batteries the caveat that those are injections or discharges.
    % Append names
    Name_vector = [Name_exist;Name_cand];

    % Table: Name vs Technology
    fuel_tech = [Thermal_exist_data.Fuel;EcoThermal_exist_data.Fuel;Renewables_exist_data.Tech;Thermal_cand_data.Fuel;EcoThermal_cand_data.Fuel;Renewables_cand_data.Tech;ESS_cand_data.Tech];

    %% Graph output
    switch graph_name_flag
        case "UCP"
            Fig_Naam = "UC-"+ opt_year +" Hourly Dispatch - "+ eval_criteria +" - sOR: " + num2str(sOR*100) + "% -ramp: " + string(v.format(ramp_capacity)) + "[MW]";
        case "GEP"
            Fig_Naam = "GEP Generation by block - "+opt_year + "-sOR: " + num2str(sOR*100) + "% -ramp: " + string(v.format(ramp_capacity)) + "[MW]";
            % For GEP we will make some markers.
            GEP_Load = GEP_Load/1000;
            x_axis_GEP_marker = [];
            for i = 1:height(GEP_Load)
                x_axis_GEP_marker = [x_axis_GEP_marker;string(v.format(GEP_Load(i)))];
            end
    end

    fig = figure("Name",Fig_Naam,'NumberTitle','off',"HandleVisibility",'on',"Visible",visibility);

    % For vertical axis readibility
    flag = 0;
    if max(sum([Output_G_exist;Output_G_cand(1:Flags.quant_cand,:)])) - min(sum(Output_G_cand(Flags.quant_cand+1:end,:))) > 20000
        Output_G_cand = Output_G_cand/1000;
        Output_G_exist = Output_G_exist/1000;
        flag = flag + 1;
    end
    
    hold on;
    title(Fig_Naam,'FontSize',10);
    % First we will graph the charging of the batteries           
    b = bar(Output_G_cand(Flags.quant_cand+1:end,:)','stacked',"FaceColor","flat");          % When GEP finished, it already gives the charging information negative, so no need to add a negative sign while making the bar graph.    
    % Color code of the bar graph that wll be below the y-axis
    for ii = 1:size(b,2)
        result = char(fuel_tech(ii));
        switch result
            case 'Coal'
                b(ii).CData = [0 0 0];                  % black
            case 'NaturalGas'
                b(ii).CData = [.5 .5 .5];               % dark gray
            case 'Bunker'
                b(ii).CData = [.25 .25 .25];            % half black half gray
            case 'Solar'
                b(ii).CData = [0.9290 0.6940 0.1250];   % opaque yellow
            case 'Wind'
                b(ii).CData = [0.3010 0.7450 0.9330];   % turquoise
            case 'Nuclear'
                b(ii).CData = [0.6350 0.0780 0.1840];   % violet
            case 'Geothermal'
                b(ii).CData = [0.8500 0.3250 0.0980];   % orange
            case 'Li'
                b(ii).CData = [.8 .8 .8];               % silver
        end
    end
    % Now we create the usual information
    b = bar([Output_G_exist;Output_G_cand(1:Flags.quant_cand,:)]','stacked',"FaceColor","flat");
    % Color code of the bar graph that will be above the y-axis
    for ii = 1:size(b,2)
        result = char(fuel_tech(ii));
        switch result
            case 'Coal'
                b(ii).CData = [0 1 1];                  % black
            case 'NaturalGas'
                b(ii).CData = [0 1 1];               % dark gray
            case 'Bunker'
                b(ii).CData = [0 1 1];            % half black half gray
            case 'Solar'
                b(ii).CData = [0 1 1];   % opaque yellow
            case 'Wind'
                b(ii).CData = [0 1 1];   % turquoise
            case 'Nuclear'
                b(ii).CData = [0 1 1];   % violet
            case 'Geothermal'
                b(ii).CData = [0 1 1];   % orange
            case 'Li'
                b(ii).CData = [0 1 1];               % silver
        end
    end
    % Figure readibility
    legend([ESS_cand_data.Name + " charging";Name_vector],"Location","bestoutside","FontSize",4);   % We will try to locate it at "best", but only because "bestoutside" cuts out the Title of the graphs
    % For vertical axis readibility
    ytickformat('%,4.0f');
    if flag == 1
        ylabel("Power [GW]");
    else
        ylabel("Power [MW]");
    end
    switch graph_name_flag
        case "UCP"
            xlabel('Time (hours)');
        case "GEP"
            xlabel('Clusters - Representative Demands [GW]');
            xticklabels(x_axis_GEP_marker)
    end
    
    ylim([min(sum(Output_G_cand(Flags.quant_cand+1:end,:))) max(sum([Output_G_exist;Output_G_cand(1:Flags.quant_cand,:)]))*1.05]);

    hold off;
    
end    

function DBB_info = call_BDD_info(Database_name)
    %% Draw all of the information from the BDD files         
    % Define existing generation data
    % Thermal
    Thermal_exist_data = readtable(Database_name,'Sheet',"ThermalExist",'ReadVariableNames',true,"UseExcel",false);
    % Variable cost According to https://www.ods.org.hn/index.php/informes/costes-marginales/2022-costomarginales/mayo22-costosmarginales:
    flag_thermal_exist = height(Thermal_exist_data);
    % EcoThermal: describing a technology that in theory is thermal, but it is not polluting (Nuclear and Geothermal)
    EcoThermal_exist_data = readtable(Database_name,'Sheet',"EcoThermalExist",'ReadVariableNames',true,"UseExcel",false);
    flag_EcoThermal_exist = flag_thermal_exist + height(EcoThermal_exist_data);
    % Renewables
    Renewables_exist_data = readtable(Database_name,'Sheet',"RenewablesExist",'ReadVariableNames',true,"UseExcel",false);
    Renewables_exist_data = sortrows(Renewables_exist_data,"Tech","ascend");
    % Suggestion: make the function to overwrite the Excel file to have it already sorted.
    flag_renewables_exist =  flag_EcoThermal_exist + height(Renewables_exist_data);
    % If this were python, we could use dictionaries and do something similar as use keys to call up specific generators
    quant_exist  = height(Thermal_exist_data) + height(EcoThermal_exist_data) + height(Renewables_exist_data);
       
    % Define candidate generation data
    % Thermal
    Thermal_cand_data = readtable(Database_name,'Sheet',"ThermalCand",'ReadVariableNames',true,"UseExcel",false);
    flag_thermal_cand = height(Thermal_cand_data);
    % EcoThermal: describing a technology that in theory is thermal, but it is not polluting (Nuclear and Geothermal)
    EcoThermal_cand_data = readtable(Database_name,'Sheet',"EcoThermalCand",'ReadVariableNames',true,"UseExcel",false);
    flag_EcoThermal_cand = flag_thermal_cand + height(EcoThermal_cand_data);
    % Renewables
    Renewables_cand_data = readtable(Database_name,'Sheet',"RenewablesCand",'ReadVariableNames',true,"UseExcel",false);
    Renewables_cand_data = sortrows(Renewables_cand_data,"Tech","ascend");
    flag_renewables_cand = flag_EcoThermal_cand + height(Renewables_cand_data);
    % Energy Storage Systems
    % Batteries
    ESS_cand_data = readtable(Database_name,'Sheet',"ESSCand",'ReadVariableNames',true,"UseExcel",false);
    flag_ESS_cand = flag_renewables_cand + height(ESS_cand_data);
    % If this were python, we could use dictionaries and do something similar as use keys to call up specific generators
    quant_cand  = height(Thermal_cand_data) + height(Renewables_cand_data) + height(EcoThermal_cand_data) + height(ESS_cand_data);
    
    % Define the fuel price variation data
    Fuel_price_data = readtable(Database_name,'Sheet',"FuelPrice",'ReadVariableNames',true,"UseExcel",false);
    RowNaams = ["2021"];
    for i = 2:30
        RowNaams(i) = [convertCharsToStrings(num2str(2021-1+i))];
    end
    Fuel_price_data.Properties.RowNames = RowNaams;

    %% Organize all flags in a single table
    Flags = table(flag_thermal_exist,flag_EcoThermal_exist,flag_renewables_exist,quant_exist,flag_thermal_cand,flag_EcoThermal_cand,flag_renewables_cand,flag_ESS_cand,quant_cand);

    %% Close all Excel files possibly still open
    fclose('all');

    DBB_info = struct('Thermal_exist_data',Thermal_exist_data,'EcoThermal_exist_data',EcoThermal_exist_data,'Renewables_exist_data',Renewables_exist_data,'Thermal_cand_data',Thermal_cand_data,'EcoThermal_cand_data',EcoThermal_cand_data,'Renewables_cand_data',Renewables_cand_data,'ESS_cand_data',ESS_cand_data,'Fuel_price_data',Fuel_price_data,'Flags',Flags);
end

function GEP_Load_tbl = calc_GEP_Load(Global_LDC_Demand,opt_LDC,fst_study_horizon,last_study_horizon)
    %Extract the headers of Global_LDC_Demand variable
    header = Global_LDC_Demand.Properties.VariableNames';
    % Store in a vector
    horizon = nan(length(header),1);
    % Convert it to number to be able to make comparisons in the following if decision maker
    for i=1:length(header)
        horizon(i) = str2double(header(i));
    end
    % Create the matrix storing the load duration and the centroids
    for j = 1:length(horizon)
        if j == 1
            GEP_Load_mat = opt_LDC(:,1); % In the first iteration we take the load duration of the LDC curve.
        else
            if fst_study_horizon > horizon(j)       % Skip years that are below the first year of the study horizon
                continue
            elseif last_study_horizon < horizon(j)  % Terminate the for loop in case we are over the last year of the study horizon
                break
            else                                    % Main body of the for loop. Appends the centroids to the GEP Load.
                GEP_Load_mat = [GEP_Load_mat,unique(Global_LDC_Demand.(num2str(horizon(j))),"stable")];
            end            
        end
    end
    
    
    
    % Convert it to Table
    % VariableNames
    header = ["LD"];
    current_year = fst_study_horizon;
    % Create the header vector
    while current_year <= last_study_horizon
        header = [header;num2str(current_year)];
        current_year = current_year + 1;
    end
    
    % Cannot convert directly to table by: GEP_Load = table(GEP_Load,'VariableNames',header);
    % So, we will create another for loop to be able to append the columns.
    
    for i = 1:length(header)
        if i == 1
            GEP_Load_tbl = table(GEP_Load_mat(:,i),'VariableNames',header(i));
        else
            GEP_Load_tbl = [GEP_Load_tbl,table(GEP_Load_mat(:,i),'VariableNames',header(i))];
        end        
    end
end

function graph_GEP(Output_G_exist,Output_G_cand,iteration,fst_study_horizon,visibility,GEP_Load,DBB_info,sOR,ramp_capacity)        
% Documentation: https://undocumentedmatlab.com/articles/export_fig
%                https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig?s_tid=mwa_osa_a
% Read ME:       https://github.com/altmany/export_fig/blob/master/README.md
% Ghostscript:   https://www.mathworks.com/matlabcentral/answers/252411-how-to-actually-use-export_fig

% Basic syntax: export_fig(filename, [handle], options...)
    Route = pwd + "\\Reports";     % Directory name

    for t = 1:size(Output_G_exist,3)    % We then create all of the figures
        opt_year = string(fst_study_horizon+t-1);
        fig = graph_EnergyProduction(Output_G_exist(:,:,t),Output_G_cand(:,:,t),opt_year,"GEP",sOR(t),"",visibility,GEP_Load.(opt_year),DBB_info,ramp_capacity(t));
        % [imageData, alpha] = export_fig:          we disregard imageData and alpha because there is no need for it.
        [~,~] = export_fig(Route + "\GEP_iteration-"+ num2str(iteration) +"_sOR-"+ num2str(max(sOR)*100)+"_ramp "+ num2str(max(ramp_capacity)) +"MW.pdf",[fig],'-append');           % In this step we export the figures to pdf. ,'-depsc'
    end

    %% Close all Excel files possibly still open
    fclose('all');      close all hidden;       close all force;
end

function graph_UCP(Output_G_exist,Output_G_cand,sOR,opt_year,eval_criteria,visibility,DBB_info,ramp_capacity)          %%FIX: This module should use export_fig(). However, for some reason due to Ghoscript internal errors, the Script is not working. I am using right now an ugly print setup.
    Route = pwd + "\\Reports";     % Directory name

    fig = graph_EnergyProduction(Output_G_exist,Output_G_cand,opt_year,"UCP",sOR,eval_criteria,visibility,[],DBB_info,ramp_capacity);     % We create a figure which displays the UCP for a given week. 
    exportgraphics(fig,Route + "\\UCP-year-"+ opt_year +".pdf","Append",true,'Resolution',900);      % https://www.mathworks.com/help/matlab/ref/exportgraphics.html#namevaluepairarguments

    %% Close all Excel files possibly still open
    fclose('all');      close all hidden;       close all force;
end

function [yearly_Installed_Capacity,legend_description] = graph_Peak_Invest(Investment,GEP_Load,flag_peakLoad_InstalledCapacity,yearly_Installed_Capacity,legend_description,renew_gen_block,blocks,DBB_info,year,sOR,ramp_capacity)
    % For readibility
    import java.text.*;     v = DecimalFormat;  % To use when printing results for easier readibility.
    % 1. First, we need to call the BDD info to make sure we are using the correct placement of all info
    Thermal_exist_data = DBB_info.Thermal_exist_data;
    EcoThermal_exist_data = DBB_info.EcoThermal_exist_data;
    Renewables_exist_data = DBB_info.Renewables_exist_data;    
    Thermal_cand_data = DBB_info.Thermal_cand_data;
    EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
    Renewables_cand_data = DBB_info.Renewables_cand_data;
    ESS_cand_data = DBB_info.ESS_cand_data;
    % 2. We now create a vector with the installed capacity of each of the technologies
    % Existing Generator Bounds
    Exist_Installed = [Thermal_exist_data.UpperLimit; EcoThermal_exist_data.UpperLimit; Renewables_exist_data.UpperLimit]; % Upper Generation Bounds Existing
    % Candidate Generator Bounds
    Cand_Installed = [Thermal_cand_data.UpperLimit; EcoThermal_cand_data.UpperLimit; Renewables_cand_data.UpperLimit;ESS_cand_data.UpperLimit]; % Upper Generation Bounds Candidate
    % 3. Multiply the Installed capacities by their respective peak demand contribution factors
    % Peak Factor coefficient factor contribution for existing generators
    % Renewables Peak Factor coefficient factor
    Exist_renewables_capacity_factor = ones(blocks,height(Renewables_exist_data));
    for i = 1:height(Renewables_exist_data)
        Exist_renewables_capacity_factor(:,i) = [renew_gen_block.(string(Renewables_exist_data.RenewableScenario{i}))];
    end
    % Take out the capacity factor of renewables for peak demand
    pPC_exist = [Thermal_exist_data.HistoricalAvailability;EcoThermal_exist_data.HistoricalAvailability;Exist_renewables_capacity_factor(1,:)']; % Peak factor of all technologies in peak demand.
    
    % Peak Factor coefficient factor contribution for candidate generators
    Cand_renewables_capacity_factor = ones(blocks,height(Renewables_cand_data));
    for i = 1:height(Renewables_cand_data)
        Cand_renewables_capacity_factor(:,i) = [renew_gen_block.(string(Renewables_cand_data.RenewableScenario{i}))];
    end
    % Take out the capacity factor of renewables for peak demand
    pPC_cand = [Thermal_cand_data.HistoricalAvailability;EcoThermal_cand_data.HistoricalAvailability;Cand_renewables_capacity_factor(1,:)';ones(height(ESS_cand_data),1)]; % Peak factor of all technologies in peak demand.        We define the factor of Batteries as able to deliver all of their installed capacity.

    % Multiply Exist_Installed.*pPC_exist and also Can_Installed.*pPC_cand
    FirmPower_Exist = Exist_Installed.*pPC_exist;       % Existing
    FirmPower_Cand  = Cand_Installed.*pPC_cand;         % Candidate

    % Multiply vector FirmPower_Cand for each of the years in Investment
    yearly_FirmPower_Cand = FirmPower_Cand'*Investment;     % this results in a vector in which each entry is the firm power for each year
    
    % 4. We will make some tags for the legend that differentiates each iteration with sOR and ramp at year it failed.
    legend_description = [legend_description; "year: "+ string(year)+ " sOR: " + string(sOR*100) + "% ramp: "+ string(v.format(ramp_capacity))];
    yearly_Installed_Capacity = [yearly_Installed_Capacity; sum(FirmPower_Exist) + yearly_FirmPower_Cand];          % This vector stores the yearly Firm Power per iteration

    % Now we decide to graph or not. In case we do not graph we will update a table called yearly Installed capacity that will store the information for all iterations
    if flag_peakLoad_InstalledCapacity == 1
        % Also, we need to extract the peak demand for each year
        PLoad = GEP_Load{1,2:end}; % We only extract peak demand
        
        % We extract the years being graphed
        years = GEP_Load(1,2:end);                  % we extract the information as table
        years = years.Properties.VariableNames;     % we now have the years stored in a cell char array
        label = years;                              % we use this helper variable called label
        years = [];
        for ii = 1:length(label)
            years = [years; str2double(cell2mat(label(ii)))];
        end
    
        fig = figure("Name","Yearly Peak Demand and Firm Power per iteration",'NumberTitle','off',"HandleVisibility",'on',"Visible",'on');     
        
        % First graph is the bar plot demand
        hold on;
        title("Peak Demand and Firm Power per iteration");          
        axis([-inf inf 0 inf]);                                 % We define the floor of the graph as zero. The rest are automatic.
        if max(max(yearly_Installed_Capacity)) >= 30000 || max(PLoad) >= 30000
            yearly_Installed_Capacity = yearly_Installed_Capacity/1000;
            PLoad = PLoad / 1000;
            ylabel('Power [GW]')
        else
            ylabel('Power [MW]')
        end
        bar(years,PLoad);                         % We make that each line are marked with an asterisk
        
        % Rest of the graphs are the line plots for each Firm Power demand
        % Plot the iterations Installed Capacity
        for jj = 1:height(yearly_Installed_Capacity)
            plot(years,yearly_Installed_Capacity(jj,:),'-*')
        end
        
        % Figure readibility
        % Change legend description to fit in less information
        divisor = ceil(length(legend_description)/10);
        if divisor >= 2
            new_legend_index = length(legend_description):-divisor:1;
            if ~isempty(find(new_legend_index,1))
                new_legend_index = [new_legend_index,1];
            end
            new_legend_index = sort(new_legend_index,'ascend');
    
            % Correct legend description to only allocate 
            for kk = 1:length(legend_description)
                if isempty(find(new_legend_index == kk))
                    legend_description(kk) = '';
                end
            end
        end
    
        % Add Legend
        legend(["Peak Demand";legend_description],"Location","bestoutside");
        xlabel('Years');
        ytickformat('%,4.0f');
        hold off
        
        % Export the fig file
        Route = pwd + "\\Reports";     % Directory name
        % [imageData, alpha] = export_fig:          we disregard imageData and alpha because there is no need for it.
        [~,~] = export_fig(Route + "\Yearly Peak Demand and Firm Power per iteration.pdf",[fig],'-append');           % In this step we export the figures to pdf.
    end
end

function ramp = ramp_calc(Investment,iterations,DBB_info)
    % Call in information to be used to calculate information
    Thermal_exist_data = DBB_info.Thermal_exist_data;
    EcoThermal_exist_data = DBB_info.EcoThermal_exist_data;   
    Thermal_cand_data = DBB_info.Thermal_cand_data;
    EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
    ESS_cand_data = DBB_info.ESS_cand_data;
    Flags = DBB_info.Flags;
    
    % Create the ramping variable for existing technologies
    exist_ramp = 0;
    % Add the ramping capacities for existing thermal and EcoThermal
    exist_ramp = exist_ramp + sum(Thermal_exist_data.RampUp) + sum(EcoThermal_exist_data.RampUp);
    
    % Calculate the cand_ramp for year evaluated
    cand_ramp = 0;  % Initialize the ramp at zero
    if iterations ~= 0  % Cand_ramp = Thermal Ramp + EcoThermal Ramp + Battery Installed Capacity. We consider that batteries have instantaneous ramps, so it is only limited by their Installed Capacity  
        cand_ramp = cand_ramp + Thermal_cand_data.RampUp'*Investment(1:Flags.flag_thermal_cand) + EcoThermal_cand_data.RampUp'*Investment(Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand) + ESS_cand_data.MaxInstall'*Investment(Flags.flag_renewables_cand+1:Flags.flag_ESS_cand);
    end
    
    % Total ramp
    ramp = exist_ramp + cand_ramp;    
end

function Renew_Cap = renew_cap_calc(Investment,renew_gen_block,iterations,DBB_info,GEP_Load,Renew_Cap,t)
    if iterations == 0
        Renew_Cap = 8760*GEP_Load.LD'*table2array(GEP_Load(:,2:end));
    else
        % Take out informaiton from DBB_info
        Thermal_exist_data = DBB_info.Thermal_exist_data;
        EcoThermal_exist_data = DBB_info.EcoThermal_exist_data;
        Renewables_exist_data = DBB_info.Renewables_exist_data;    
        Thermal_cand_data = DBB_info.Thermal_cand_data;
        EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
        Renewables_cand_data = DBB_info.Renewables_cand_data;
        ESS_cand_data = DBB_info.ESS_cand_data;
        Flags = DBB_info.Flags;
        % Existing Generator Bounds
        Exist_Upper = [Thermal_exist_data.UpperLimit; EcoThermal_exist_data.UpperLimit; Renewables_exist_data.UpperLimit]; % Upper Generation Bounds Existing
        % Candidate Generator Bounds
        Cand_Upper = [Thermal_cand_data.UpperLimit; EcoThermal_cand_data.UpperLimit; Renewables_cand_data.UpperLimit; ESS_cand_data.UpperLimit]; % Upper Generation Bounds Candidate

        
        Measure_Cap = 0;
        % Existing Renewables
        counter = 1;
        for i = Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist
            Measure_Cap = Measure_Cap + sum(Exist_Upper(i).*renew_gen_block.(string(Renewables_exist_data.RenewableScenario{counter})));
            counter = counter +1;
        end
        % Candidate Renewables
        counter = 1;
        for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
            Measure_Cap = Measure_Cap + sum(Investment(i)*Cand_Upper(i).*renew_gen_block.(string(Renewables_cand_data.RenewableScenario{counter})));
            counter = counter + 1;
        end
        % Measurement of Renew Cap
        Renew_Cap(t) = Measure_Cap;
    end
end

%% Content