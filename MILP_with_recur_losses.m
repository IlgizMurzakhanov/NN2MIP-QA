clear; clear all; clc; close all;
% MILP implementation with optional epsilon reduction
% eps = 0 3650.21 $/MW
% eps = 10 3667.237 -> stable
% Pslack explicitly (sdpvar)
% constants
epsR = 10;%4.125;
warning('off')
verbose = false;
log = false;%false;
logfile = strcat('MILP_PS_',datestr(now,'dd_mm_yyTHH_MM_SS'),'.log');%'out_1203.log';

% mpc = util.case14_wind;     % base mpc case
% mpc = util.case14;     % base mpc case
mpc = util.case14_wind_corr;     % base mpc case
%% relax q limits
% QMAX = 4;
% QMIN = 5;
% mpc.gen(:,QMAX) = mpc.gen(:,QMAX)+0.25*(mpc.gen(:,QMAX)-mpc.gen(:,QMIN));
% mpc.gen(:,QMIN) = mpc.gen(:,QMIN)-0.25*(mpc.gen(:,QMAX)-mpc.gen(:,QMIN));
% 
% sTol = 1.3;
% mpc.branch(:,6) = mpc.branch(:,6)*sTol; 
% mpc.branch(:,7) = mpc.branch(:,7)*sTol;
% mpc.branch(:,8) = mpc.branch(:,7)*sTol;
%%
mpcOpt = mpoption;
mpcOpt.verbose = 0;
mpcOpt.out.all = 0;
qOpt = 0;
mpcOpt.pf.enforce_q_lims = qOpt;

%path = strcat(pwd,filesep,'neural_network',filesep,'balanced_set',filesep);
%path = '~/Documents/powersys/neural_network/balanced_set/'; % NN formulation
%path = '~/Documents/powersys/power_system_database_generation/neural_network/test_set2/'; % NN formulation
% path = [pwd,filesep,'neural_network',filesep,'qLimitsNotEnforcedQS3_2relaxed',filesep]; % dynamic path to q lims not enforced, relaxed
path = 'C:\Users\������\Documents\���������_������\DTU\Research\Collaboration with Andreas\neural_network\qLimitsNotEnforcedQS3relaxed\';
path = 'C:\Users\������\Documents\���������_������\DTU\Research\Collaboration with Andreas\neural_network\qLimitsNotEnforced\';
path = 'C:\Users\������\Documents\���������_������\DTU\Research\Collaboration with Andreas\neural_network\14_bus_system_June2021\14_bus_PG_variables\';
% path = 'C:\Users\������\Documents\���������_������\DTU\Research\Collaboration with Andreas\neural_network\balanced_set\'; %did not end in 240 iterations
HIDDEN_LAYER = 50;          % MILP formulation
RELU_LAYERS = 3;            % MILP formulation
DELTA_P = 10;               % losses

%% output and log
if verbose
   solverOutput = {'verbose',1}; %#ok
else
   solverOutput = {'verbose',0};
end

if log
    logFile = strcat(pwd,filesep,'log',filesep,logfile); %#ok
    fprintf('Logging all output to:\n%s\n',logFile)
    fprintf('Open <a href="matlab: opentoline(''%s'',1)">file</a>!\n',logFile)
    diary on
    diary(logFile)
    fprintf(['\n\n[',datestr(now),'] - ',mfilename,'\n'])
else
    diary off
end

%% -----  NN data  -----
[W_input,W,W_output,bias, NN_input, NN_output] = milp.get_nn_data(path);

%% -----  inverval arithmetic  -----

interval_arithmetic = true;

x_0_up = ones(HIDDEN_LAYER,1,RELU_LAYERS)*(1000);
x_0_lp = ones(HIDDEN_LAYER,1,RELU_LAYERS)*(-1000);

if (interval_arithmetic == true)
    u_init = ones(size(W_input, 2),1);
    l_init = zeros(size(W_input,2),1);
    x_0_up(:,1,1) = max(W_input,0)*u_init+min(W_input,0)*l_init+bias{1};
    x_0_lp(:,1,1) = min(W_input,0)*u_init+max(W_input,0)*l_init+bias{1};
    for j = 1:RELU_LAYERS-1
        x_0_up(:,1,j+1) = max(W{j},0)*max(x_0_up(:,1,j),0)+min(W{j},0)*max(x_0_lp(:,1,j),0)+bias{j+1};
        x_0_lp(:,1,j+1) = min(W{j},0)*max(x_0_up(:,1,j),0)+max(W{j},0)*max(x_0_lp(:,1,j),0)+bias{j+1};
    end
end


%% -----  milp formulation  -----
% construct otpimization problem of neural network

size_input = size(NN_input,2);
LP_relax = false;

u_NN = sdpvar(size_input,1); %scaled between 0 - 1 w/ generator range
Pslack = sdpvar(1,1); % PU


if LP_relax == true     % integer relaxation
    ReLU_0 = sdpvar(HIDDEN_LAYER,1,RELU_LAYERS);
else
    ReLU_0 = binvar(HIDDEN_LAYER,1,RELU_LAYERS);
end

x_0 = sdpvar(HIDDEN_LAYER,1,RELU_LAYERS);
x_0_ReLU = sdpvar(HIDDEN_LAYER,1,RELU_LAYERS);
y = sdpvar(2,1);

% ----- objective function -----
genRange = [60 60 25 25]';
optvar2pu = genRange./mpc.baseMVA;
P_d = sum(mpc.bus(:,3));

obj = mpc.gencost(2:end,6)'*(u_NN.*genRange)+mpc.gencost(1,6)*Pslack*mpc.baseMVA;

i_loss = 0;

mpcOrigin = mpc;
err = 1;
  %% ----- loss approximation -----
[l,pGeneration] = milp.getLossSensitivity(mpc,mpcOpt,DELTA_P,solverOutput{:});
while err >= 0.01   

    % ----- constraint -----

    [PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
    constr = [];
    constr = [constr;...
        0<=u_NN<=1]:'input restriction';%#ok % input restrictions
    % ----- Slack bus -----
    constr = [constr;...
        (0 <= Pslack <= 6.15):'slack bus']; %#ok

    P_demand_i = sum(mpc.bus(:,PD))/mpc.baseMVA;% change to PD - sum(PG) < PslackMax

    constr = [constr;...
        (Pslack + sum(u_NN.*optvar2pu) - P_demand_i - l(1)...
        - l(2)*(u_NN(1)*optvar2pu(1) - pGeneration(2)) ...
        - l(3)*(u_NN(2)*optvar2pu(2) - pGeneration(3))...
        - l(4)*(u_NN(3)*optvar2pu(3) - pGeneration(4))...
        - l(5)*(u_NN(4)*optvar2pu(4) - pGeneration(5)) == 0):'lossy power balance'];%#ok
    constr = [constr; ...
        (x_0(:,:,1) == W_input*u_NN + bias{1}):'input ReLU']; %#ok

    for i = 1:RELU_LAYERS
        constr = [constr; ...
            (x_0_ReLU(:,:,i) <= x_0(:,:,i) - x_0_lp(:,:,i).*(1-ReLU_0(:,:,i))):'hidden ReLU 1';...
            (x_0_ReLU(:,:,i) >= x_0(:,:,i)):'hidden ReLU 2';...
            (x_0_ReLU(:,:,i) <= x_0_up(:,:,i).*ReLU_0(:,:,i)):'hidden ReLU 3';...
            (x_0_ReLU(:,:,i) >= 0):'hidden ReLU 4'];%#ok
    end
    for i = 1:RELU_LAYERS-1
        constr = [constr; ...
            (x_0(:,:,i+1) == W{i}*x_0_ReLU(:,:,i) + bias{i+1}):'hidden ReLU eq'];%#ok
    end

    % ----- integer relaxation -----
    if LP_relax == true
        % integer relaxation
        constr = [constr; ...
            (0<= ReLU_0 <=1):'int relax' ];%#ok
    end

    %% ----- output layer -----
    constr = [constr; ...
        (y == W_output * x_0_ReLU(:,:,end) + bias{end}):'output ReLU'];%#ok % I am not sure about this..
    if(epsR) 
        constr = [constr; ...
            (y(2) >= y(1)+epsR):'classification'];%#ok;]; % classification as unsafe
    else
        constr = [constr; ...
            (y(2) >= y(1)):'classification'];%#ok;]; % classification as unsafe
    end
    %% ----- optimization -----
    % set Gurobi as solver
    options = sdpsettings('solver','GUROBI',solverOutput{:});
    % solve
    diagnostics = optimize(constr,obj,options);

    % ----- results -----
    fprintf('\nSolution found at: %2.3f\t\n',value(obj))
    setpoints = [value(Pslack)*mpc.baseMVA; value(u_NN).*genRange];
    milpRes.setpoints = setpoints;
    milpRes.mpc = mpc;
    fprintf('Optimal solution found at: \n\t  MW \t   $/MW\n')
    for i = 1:5
        fprintf('\t%2.3f\t | %2.1f\n',setpoints(i), mpc.gencost(i,6));
    end
    fprintf('\n')
    % setpoints = value(u_NN).*genRange;

    [err,~,pfSol] = getLosses(milpRes, mpcOpt);
     fprintf('Loss estimation error: [%d] - %2.3f\n',i_loss,err);

    [l,pGeneration] = iterateLosses(mpc, mpcOpt, pfSol,DELTA_P);
    mpc = pfSol;
    i_loss = i_loss + 1;
end



%% ----- check ACOPF solution -----
setpoints = value(u_NN).*genRange;
acopfRes = analyzeResults(mpcOrigin,mpcOpt,obj,setpoints,genRange,NN_input, NN_output);
fprintf('Loss iterations: %d\n',i_loss)
fprintf('Loss estimation error: %2.3f\n',err);
if epsR
    fprintf('Epsilon reduction applied <strong>eps=%3.2f</strong>\n',epsR) %#ok
end
diary off
if log
    fprintf('Open <a href="matlab: opentoline(''%s'',1)">log file</a>!\n',logFile) %#ok
end











%% ----- functions -----
function [l,pGeneration] = iterateLosses(mpc,mpopt,pfResults,DELTA_P)
PG = 2;
PD = 3;
    pGeneration = pfResults.gen(:, PG);
    genlist = util.GetGenList(mpc); % 0 for slack 1 for rest
    NG = length(genlist);
    l = nan(NG,1);

    %% 1 st step -> calculate Pi_DC
    l(1) = sum(pGeneration) - sum(mpc.bus(:, PD)); % PD = 3 real power demand

    %% 2 nd step -> ACOPF to calculate P_L^0, \gamma
    for i = 2:NG
        mpcTemp = mpc;
        mpcTemp.gen(:, PG) = pGeneration;
        mpcTemp.gen(i, PG) = pGeneration(i)+ DELTA_P;

        pf = runpf(mpcTemp, mpopt);
        loss = sum(pf.gen(:, PG)) - sum(mpc.bus(:, PD));
        l(i) = (loss - l(1))/DELTA_P; % don't divide this one with baseMVA, is delta P is given in MW
    end
    pGeneration = pGeneration./mpc.baseMVA;
    l(1) = l(1)/mpc.baseMVA;
end

function [err, varargout] = getLosses(milpRes,mpOpt)
% get loss comparision of MILP and PF results
% return loss estimation error
% optionally return loss values

    PD = 3;
    PG = 2;
    PMAX = 9;
    PMIN = 10;
    TOL = 1e-2;
    % References: https://matpower.org/docs/ref/matpower7.0/lib/get_losses.html
    %% Prepare AC PF to run with MILP setpoints
    milpRes.mpc.gen(:,PMAX) = milpRes.setpoints + TOL;
    milpRes.mpc.gen(:,PMIN) = milpRes.setpoints - TOL;
    milpRes.mpc.gen(:,PG) = milpRes.setpoints;
    
    pfRes = runpf(milpRes.mpc,mpOpt);

    loss.ac_real = sum(real(get_losses(pfRes)));
    loss.milp_real = sum(milpRes.setpoints) - sum(milpRes.mpc.bus(:,PD));
    err = abs(loss.ac_real - loss.milp_real)/loss.ac_real;

    if nargout >= 2
        varargout{1} = loss;
    end
    if nargout >= 3
        varargout{2} = pfRes;
    end
end

function varargout = analyzeResults(mpc,mpOption,obj,setpoints,genRange,NN_input, NN_output)
chr = fprintf('|  MILP objective value:     <strong>%3.2f $/MW </strong> |\n',value(obj));
chr = chr -17;
fprintf(strcat('|', repmat('-',1,chr-3),'|\n'))

genList = util.GetGenList(mpc);
[classification,classDetails,dampingR,pfSuccess] = util.classifyCase(mpc,setpoints,genList,mpOption);

pretty_print(pfSuccess,classification,chr)
fprintf('MILP details:\n')
fprintf('Damping ratio:\t%3.2f\n',dampingR)
fprintf('   | P  |  Q |  Vm  |  S  | SSS | N_contingency\n')
disp(classDetails)

% ----- compare ACOPF price -----
% In order to compare the MILP results w/ acopf result, the obj has to be
% scaled to $/hr. The first step to do so is to include the Slack bus
% price.
[results,~] = runopf(mpc,mpOption);
fprintf('|  ACOPF objective value:   <strong> %3.2f $/MW </strong> |\n',results.f);
fprintf(strcat('|', repmat('-',1,chr-3),'|\n'))

setpointsACOPF = results.var.val.Pg(2:end) * results.baseMVA;
[classification,classDetails,dampingR,pfSuccess] = util.classifyCase(mpc,setpointsACOPF,genList,mpOption);

pretty_print(pfSuccess,classification,chr)
fprintf('ACOPF details:\n')
fprintf('Damping ratio:\t%3.2f\n',dampingR)
fprintf('   | P  |  Q |  Vm  |  S  | SSS | N_contingency\n')
disp(classDetails)
%% ----- distance from feasibility -----
compare_feas_dist(setpoints, setpointsACOPF,genRange,NN_input, NN_output,chr)
if nargout ==1 
   varargout{1} = results;
end
end

function pretty_print(pfSuccess,classification,chr)
fprintf('|\tPowerFlow convergence:\t\t\t\t%d  |\n',pfSuccess)
fprintf('|\tACOPF constraints:\t\t\t\t\t%d  |\n',classification(1))
fprintf('|\tDamping Ratio constraints:\t\t\t%d  |\n',classification(2))
fprintf(strcat('|', repmat('=',1,chr-3),'|\n'))
end

function compare_feas_dist(setpoints, setpointsACOPF,genRange,NN_input, NN_output,chr)
% Calculate distance from nearest feasible point. Nearest feasible point is
% calculated from NN_input and NN_output, where input is an array of
% setpoints and output is the classification correspoinding to the
% setopints.
fprintf('\n Distances\n')
fprintf(strcat('|', repmat('-',1,chr-3),'|\n'))
acopf_milp = sum(abs(setpoints-setpointsACOPF));
fprintf('|\tACOPF - MILP:\t\t<strong>%3.2f</strong>\t\t\t   |\n',acopf_milp)

allStablePoints = NN_input(logical(NN_output(:,1)),:).*genRange';
[minVal, minIdx] = min(sum(abs(allStablePoints - setpoints'),2));
fprintf('|Distance from closest feasible point:\t   |\n|\t<strong>MILP\t-> \t\t%2.3f</strong>\t\t\t\t   |\n',minVal)
fprintf('|DB:\t%2.3f  %2.3f  %2.3f  %2.3f\t   |\n',allStablePoints(minIdx,:));
fprintf('|MILP:\t%2.3f  %2.3f  %2.3f  %2.3f\t   |\n',setpoints)

[minVal, minIdx] = min(sum(abs(allStablePoints - setpointsACOPF'),2));
fprintf('|\t<strong>ACOPF\t-> \t\t%2.3f\t</strong>  \t\t\t   |\n',minVal)
fprintf('|DB:\t%2.3f  %2.3f  %2.3f  %2.3f\t   |\n',allStablePoints(minIdx,:));
fprintf('|ACOPF:\t%2.3f  %2.3f  %2.3f  %2.3f\t   |\n',setpointsACOPF)
fprintf(strcat('|', repmat('=',1,chr-3),'|\n\n'))
end
