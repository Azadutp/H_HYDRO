
clear all 
clc
close all

%% 1 - Load data
% Two variables are loaded from data set:
%   inFlow - flow rate into the reservoir measured hourly (CFS)
%   price - price of electricity measured hourly ($/MWh)
load('FlowAndPriceData.mat');



%% 2 - Define Starting Point and Constants
% We wish to find the vector x of optimal decision variables. 
% For this problem:
% x = [turbineFlow(1) ... turbineFlow(N) spillFlow(1) ... spillFlow(N)]

% Our initial starting point is turbineFlow = inFlow and spillFlow = 0
x0 = [inFlow; zeros(size(inFlow))];

stor0 = 90000; % initial vol. of water stored in the reservoir (Acre-Feet)

k1 = 0.00003; % K-factor coefficient
k2 = 9; % K-factor offset

MW2kW = 1000; % MW to kW
C2A = 1.98347/24; % Convert from CFS to AF/HR

C = [C2A k1 k2 MW2kW]; % Vector of constants

N = length(inFlow);

%% 3 - Define Linear Inequality Constraints and Bounds
DefineConstraints;


%% 4 - Objective function
% Use symbolic math for creating the objective function and then convert
% into a MATLAB function that can be evaluated.
%% 4.1 - Derive Objective Function using Symbolic Math
% First define the decision variables as X, which is a vector consisting
% of [TurbineFlow, SpillFlow] and all the parameters as symbolic objects
X = sym('x',[2*N,1]);   % turbine flow and spill flow
F = sym('F',[N,1]);     % flow into the reservoir
P = sym('P',[N,1]);     % price
s0 = sym('s0');         % initial storage
c = sym({'c1','c2','c3','c4'}.','real');   % constants

%%
% Total out flow equation
TotFlow = X(1:N)+X(N+1:end);

%%
% Storage Equations
S = cell(N,1);
S{1} = s0 + c(1)*(F(1)-TotFlow(1));

for ii = 2:N
    S{ii} = S{ii-1} + c(1)*(F(ii)-TotFlow(ii));
end

%%
% K factor equation
k = c(2)*([s0; S(1:end-1)]+ S)/2+c(3);

% MWh equation
MWh = k.*X(1:N)/c(4);

%%
% Revenue equation
R = P.*MWh;

%%
% Total Revenue (Minus sign changes from maximization to minimization)
totR = -sum(R);

%% 4.2 - Create Function File that can be Evaluated in MATLAB
% Create a function that when evaluated returns the value at a point

% Perform substitutions so that X is the only remaining variable 
totR = subs(totR,[c;s0;P;F],[C';stor0;price;inFlow]);

% Generate a MATLAB file from the equation
matlabFunction(totR,'vars',{X},'file','objFcn');


%%

SearchAgents_no=50; % Number of search agents

Function_name='F1'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)

Max_iteration=500; % Maximum numbef of iterations

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

[Best_score,Best_pos,GWO_cg_curve]=GWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

[Best_score_CS,Best_pos_CS,GWOCS_cg_curve]=GWO_CS(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);


figure %('Position',[500 500 660 290])
%Draw search space
subplot(1,2,1);
func_plot(Function_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])

%Draw objective space
subplot(1,2,2);
semilogy(GWO_cg_curve,'Color','r')
hold on
semilogy(GWOCS_cg_curve,'Color','k')
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');

axis tight
grid on
box on
legend('GWO','GWOCS')

display(['The best solution obtained by GWO is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective funciton found by GWO is : ', num2str(Best_score)]);
display(['The best optimal value of the objective funciton found by GWOCS is : ', num2str(Best_score_CS)]);

        



