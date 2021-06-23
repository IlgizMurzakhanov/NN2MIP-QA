clear; clear all; clc; close all;

% Solve N MILP's change gen prices (spaceVec) and eps values (epsArray),
% Pslack explicitly (sdpvar)
% Dependencies: install GUROBI and YALMIP and IPOPT
%
% ITERATION ON EPS ARRAY is neccessary in order to get information about
% the share of stable points.
% util.epsChangeSummary - evaluate results
tic
warning('off')
log = true;%false;
epsArray = 0;
fileNamePre = strcat('MILP_PRICE_CHANGES_QS3_',datestr(now,'dd_mm_yyTHH_MM_SS'));
logfile = strcat(fileNamePre,'.log');%'out_1203.log';
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
% constants
mpc = util.case14_wind;     % base mpc case

%% modify mpc
% ====== relax q limits
QMAX = 4;
QMIN = 5;
mpc.gen(:,QMAX) = mpc.gen(:,QMAX)+0.25*(mpc.gen(:,QMAX)-mpc.gen(:,QMIN));
mpc.gen(:,QMIN) = mpc.gen(:,QMIN)-0.25*(mpc.gen(:,QMAX)-mpc.gen(:,QMIN));
% ===== Relax S limits
sTol = 1.3;
mpc.branch(:,6) = mpc.branch(:,6)*sTol; 
mpc.branch(:,7) = mpc.branch(:,7)*sTol;
mpc.branch(:,8) = mpc.branch(:,7)*sTol;

%% mpopt
mpcOpt = mpoption;
mpcOpt.verbose = 0;
mpcOpt.out.all = 0;

qOpt = 0;
mpcOpt.pf.enforce_q_lims = qOpt;
solverOutput = {
    'verbose',0};

% path = [pwd,filesep,'neural_network',filesep,'qLimitsNotEnforcedQS3_2relaxed',filesep];
path = 'C:\Users\Ильгиз\Documents\Документы_Ильгиз\DTU\Research\Collaboration with Andreas\neural_network\qLimitsNotEnforcedQS3relaxed\';

HIDDEN_LAYER = 50;          % MILP formulation
RELU_LAYERS = 3;            % MILP formulation
DELTA_P = 10;               % losses

[W_input,W,W_output,bias, NN_input, NN_output] = get_nn_data(path);

interval_arithmetic = true;

% upper bound on x_0 (Here we will need to use some bound tightening)
x_0_up = ones(HIDDEN_LAYER,1,RELU_LAYERS)*(1000);
% lower bound on x_0 (Here we will need to use some bound tightening)
x_0_lp = ones(HIDDEN_LAYER,1,RELU_LAYERS)*(-1000);

if (interval_arithmetic == true)
    % use interval arithmetic to compute tighter bounds
    % initial input bounds
    % Big M bound the solution
    % reformulation of max function
    
    u_init = ones(size(W_input, 2),1);
    l_init = zeros(size(W_input,2),1);
    x_0_up(:,1,1) = max(W_input,0)*u_init+min(W_input,0)*l_init+bias{1};
    x_0_lp(:,1,1) = min(W_input,0)*u_init+max(W_input,0)*l_init+bias{1};
    for j = 1:RELU_LAYERS-1
        x_0_up(:,1,j+1) = max(W{j},0)*max(x_0_up(:,1,j),0)+min(W{j},0)*max(x_0_lp(:,1,j),0)+bias{j+1};
        x_0_lp(:,1,j+1) = min(W{j},0)*max(x_0_up(:,1,j),0)+max(W{j},0)*max(x_0_lp(:,1,j),0)+bias{j+1};
    end
end

% construct otpimization problem of neural network
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

size_input = size(NN_input,2);
Pslack = sdpvar(1,1); % PU

LP_relax = false;

u_NN = sdpvar(size_input,1);
if LP_relax == true     % integer relaxation
    ReLU_0 = sdpvar(HIDDEN_LAYER,1,RELU_LAYERS);
else
    ReLU_0 = binvar(HIDDEN_LAYER,1,RELU_LAYERS);
end

x_0 = sdpvar(HIDDEN_LAYER,1,RELU_LAYERS);
x_0_ReLU = sdpvar(HIDDEN_LAYER,1,RELU_LAYERS);
y = sdpvar(2,1);

% OBJECTIVE FUNCTION
% Get cost functions from mpc
% mpc = util.case14_wind;

genRange = [60 60 25 25]';
optvar2pu = genRange./mpc.baseMVA;
P_demand_i = sum(mpc.bus(:,PD))/mpc.baseMVA;

%% ----- set up iterators -----
% nPoints = size(spaceGrid,1) * size(epsArray,2)

n = 1;
pMin = 0; % $/MW
pMax = 0; % $/MW
% n = 5;
% pMin = 0; % $/MW
% pMax = 5; % $/MW

pVec = linspace(pMin,pMax,n);
% ------ to change only 1 generators prices (eg pslack)
% <<<<<< comment next 2
% [pG1,pG2,pG3] = ndgrid(pVec,pVec,pVec); % meshgrid % 
% spaceGrid = [pG1(:) pG2(:) pG3(:)];
% using lhsdesign
spaceGrid = lhsdesign(length(pVec)^3,3);
spaceGrid = (spaceGrid.*(pMax-pMin))+pMin;

% <<<<<< insert next
%spaceGrid = [pVec', repmat(mpc.gencost(2,6),[n,1]),repmat(mpc.gencost(3,6),[n,1])];

nEps = size(epsArray,2);
nCostGrid = size(spaceGrid,1);
nPoints = nCostGrid * nEps;
pslack= nan(nPoints, 1);
price = nan(nPoints, 1); % price (for coloring)
stability = nan(nPoints, 1); % see if the solution did not violate constraints
dampingRatio = stability;
gsetpoints = nan(nPoints, 4);
resultStruct = cell(nPoints,3); %milp,acopf, n-1 acopf
% FIX generator 6 and 8 cost to 0

spaceGrid = [spaceGrid zeros(nCostGrid,2) ];
tic
mpcOrigin = mpc;
%% ----- START ITERATION -----
fprintf('progress:  ')
chr = fprintf('%2.2f %%\n',0);
toc
nk = 1;
for k = 1:nEps:nPoints
    obj = spaceGrid(nk,2:end)*(u_NN.*genRange)+spaceGrid(nk,1)*Pslack*mpc.baseMVA;
    % ----- loss approximation -----
    
    disp('Initial loss and pGeneration')
    [l,pGeneration] = milp.getLossSensitivity(mpc,mpcOpt,DELTA_P,solverOutput{:})
    
    for j = 0:nEps-1
        i_loss = 0;
        mpc = mpcOrigin;
        err = 1;
        %% LOSS ESTIMATION ITERATION =========================================
        while err >= 0.01  
            % ----- constraints -----
            constr = [];
            constr = [constr;...
                0<=u_NN<=1]:'input restriction';%#ok
            constr = [constr;...
                (0 <= Pslack <= 6.15):'slack bus']; %#ok
            % ----- losses -----
            constr = [constr;...
            (Pslack + sum(u_NN.*optvar2pu) - P_demand_i - l(1)...
            - l(2)*(u_NN(1)*optvar2pu(1) - pGeneration(2)) ...
            - l(3)*(u_NN(2)*optvar2pu(2) - pGeneration(3))...
            - l(4)*(u_NN(3)*optvar2pu(3) - pGeneration(4))...
            - l(5)*(u_NN(4)*optvar2pu(4) - pGeneration(5)) == 0):'lossy power balance'];%#ok
            % ----- input layer -----
            constr = [constr; ...
                (x_0(:,:,1) == W_input*u_NN + bias{1}):'input ReLU']; %#ok
            % ----- hidden layers -----
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
            % ----- output layer -----
            constr = [constr; ...
                (y == W_output * x_0_ReLU(:,:,end) + bias{end}):'output ReLU'];%#ok % I am not sure about this..
            eps = epsArray(j+1);
            constr = [constr; ...
                     (y(2) >= y(1)+eps):'classification'];%#ok;]; % with epsilon control

                    %(y(2) >= y(1)):'classification'];%#ok;]; % classification as unsafe
            % ----- optimization -----
            options = sdpsettings('solver','gurobi',solverOutput{:});
            diagnostics = optimize(constr,obj,options)
            milpRes.setpoints  = [value(Pslack)*mpc.baseMVA; value(u_NN).*genRange];
            disp('milpRes.setpoints')
            milpRes.setpoints
            
            milpRes.mpc = mpc;

            [err,~,pfSol] = getLosses(milpRes, mpcOpt);
            fprintf('Loss estimation error: [%d] - %2.3f\n',i_loss,err);

            [l,pGeneration] = iterateLosses(mpc, mpcOpt, pfSol,DELTA_P)
            
            mpc = pfSol;
            i_loss = i_loss + 1;
            toc
        end
        % ----- results -----
        setpoints = value(u_NN).*genRange;
        % ----- check solution -----
        [milpRes,acopfRes,N1Res] = util.compareResults(mpcOrigin,mpcOpt,obj,setpoints,spaceGrid(nk,:),'silent');
        pslack(k+j) = value(Pslack);
        price(k+j) = value(obj);
        stability(k+j) = milpRes.stability;
        dampingRatio(k+j) = milpRes.dr; %keep this for backwards compatibility of plots
        milpRes.time = diagnostics.solvertime;
        milpRes.diag = diagnostics;
        gsetpoints(k+j,:) = value(u_NN);%setpoints(:);
        resultStruct{k+j,1} = milpRes;
        resultStruct{k+j,2} = acopfRes;
        resultStruct{k+j,3} = N1Res;
        
        fprintf('%2.2f %%\n',(k+j)/nPoints*100)
    end
    nk = nk +1;
end

toc
warning('on')
fileName = strcat(fileNamePre,'.mat');
save(fileName,'pslack','spaceGrid','price','stability','gsetpoints','dampingRatio','resultStruct','epsArray')
%%%% saved variables
% pslack        - value of the slack generator real power output
% spaceGrid     - the grid of cost functions (each row defines a cost function)
% price         - the price of the MILP solution
% stability     - 1 if MILP satisfies all criteria
% gsetpoints    - generator setpoints returned by MILP
% damping ratio - damping ratio
%
diary off
if log
    fprintf('Open <a href="matlab: opentoline(''%s'',1)">log file</a>!\n',logFile) %#ok
end

%% ===== END =====






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
    disp('Power flow execution: success - 1, failure - 0')
    pfRes.success

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
% ==============================

function [W_input,W,W_output,bias, NN_input, NN_output] = get_nn_data(path)

path2wb = strcat(path,'output/');
% define the neural network size

% load data from Python
NN_input = csvread(strcat(path,'NN_input.csv'));
safe_log_ = csvread(strcat(path,'NN_output.csv'));      % manipulate out size
safe_log = safe_log_(:,1);
NN_output=safe_log_(:,1);

% we have to invert the weight matrices
W_input = csvread(strcat(path2wb ,'W0_p.csv')).';
W_output = csvread(strcat(path2wb,'W3_p.csv')).';
W{1} = csvread(strcat(path2wb ,'W1_p.csv')).';
W{2} = csvread(strcat(path2wb ,'W2_p.csv')).';

bias{1} = csvread(strcat(path2wb ,'b0_p.csv'));
bias{2} = csvread(strcat(path2wb ,'b1_p.csv'));
bias{3} = csvread(strcat(path2wb ,'b2_p.csv'));
bias{4} = csvread(strcat(path2wb ,'b3_p.csv'));

end



function compare_feas_dist(setpoints, setpointsACOPF,genRange,NN_input, NN_output,chr)
% Calculate distance from nearest feasible point. Nearest feasible point is
% calculated from NN_input and NN_output, where input is an array of
% setpoints and output is the classification correspoinding to the
% setopints.
fprintf('\n Distances\n')
fprintf(strcat('|', repmat('-',1,chr-3),'|\n'))
acopf_milp = sum(abs(setpoints-setpointsACOPF));
fprintf('|\tACOPF - MILP:\t\t<strong>%3.2f</strong>\t   |\n',acopf_milp)

allStablePoints = NN_input(logical(NN_output(:,1)),:).*genRange';
[minVal, minIdx] = min(sum(abs(allStablePoints - setpoints'),2));
fprintf('|Distance from closest feasible point:\t   |\n|\t<strong>MILP\t-> \t\t%2.3f</strong>\t   |\n',minVal)
fprintf('|DB:\t%2.3f  %2.3f  %2.3f  %2.3f\t   |\n',allStablePoints(minIdx,:));
fprintf('|MILP:\t%2.3f  %2.3f  %2.3f  %2.3f\t   |\n',setpoints)

[minVal, minIdx] = min(sum(abs(allStablePoints - setpointsACOPF'),2));
fprintf('|\t<strong>ACOPF\t-> \t\t%2.3f\t</strong>   |\n',minVal)
fprintf('|DB:\t%2.3f  %2.3f  %2.3f  %2.3f\t   |\n',allStablePoints(minIdx,:));
fprintf('|ACOPF:\t%2.3f  %2.3f  %2.3f  %2.3f\t   |\n',setpointsACOPF)
fprintf(strcat('|', repmat('=',1,chr-3),'|\n\n'))
end
