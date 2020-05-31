%___________________________________________________________________%
%  Grey Wolf Optimizer (GWO) source codes version 1.0               %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper: S. Mirjalili, S. M. Mirjalili, A. Lewis             %
%               Grey Wolf Optimizer, Advances in Engineering        %
%               Software , in press,                                %
%               DOI: 10.1016/j.advengsoft.2013.12.007               %
%                                                                   %
%___________________________________________________________________%

% This function containts full information and implementations of the benchmark 
% functions in Table 1, Table 2, and Table 3 in the paper

% lb is the lower bound: lb=[lb_1,lb_2,...,lb_d]
% up is the uppper bound: ub=[ub_1,ub_2,...,ub_d]
% dim is the number of variables (dimension of the problem)

function [lb,ub,dim,fobj] = Get_Functions_details(F)

        %fobj = @F1;
        lb=[-100 ];
        ub=[100 ];
        dim=100;
        
  


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
%matlabFunction(totR,'vars',{X},'file','objFcn');
%%

fobj=totR;

end
% F1

% function o = F1(x)
% 
% %o=sum(x1.^2);
% end

