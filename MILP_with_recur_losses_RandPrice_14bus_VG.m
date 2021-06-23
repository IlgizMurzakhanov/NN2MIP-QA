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
epsArray = 0:10;
fileNamePre = strcat('MILP_PRICE_CHANGES_QS3_',datestr(now,'dd_mm_yyTHH_MM_SS'));
logfile = strcat(fileNamePre,'.log');%'out_1203.log';
if log
    logFile = strcat(pwd,filesep,'log',filesep,logfile); %#ok
    fprintf('Logging all output to:\n%s\n',logFile)
    fprintf('Open <a href="matlab: opentoline(''%s'',1)">file</a>!\n',logFile)
    diary on % WAS SWITCHED ON INITIALLY
    diary(logFile)
    fprintf(['\n\n[',datestr(now),'] - ',mfilename,'\n'])
else
    diary off
end
% constants
%%%%%%%%%%%%%%%%%%%%%% SET THE CASE FILE: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mpc = util.case14_wind;     % base mpc case
% AllGenInd = mpc.gen(:,1); % 
% NGen = size(AllGenInd,1);

%% modify mpc (for 68 we will use only one created database)
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
% path = 'C:\Users\??????\Documents\?????????_??????\DTU\Research\Collaboration with Andreas\neural_network\68_bus_system\';

%%%%%%%%%%%%%%%%%%%%%% SET THE PATH: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path = 'C:\Users\Ильгиз\Documents\Документы_Ильгиз\DTU\Research\Collaboration with Andreas\neural_network\14_bus_system_PG_VG\';

HIDDEN_LAYER = 100;          % MILP formulation
RELU_LAYERS = 3;            % MILP formulation
DELTA_P = 10;               % losses % HARD-CODED: see paper
DELTA_V = 0.05;

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

%% mpc case data preprocessing
% genRange
GenBusInd = mpc.bus(mpc.bus(:,2) == 2); % indices of gen-s in bus table (EXCLUDING slack node)
[~,GenGenInd] = intersect(mpc.gen(:,1),GenBusInd); % .. in gen table ..
genRange = (mpc.gen(GenGenInd,9) - mpc.gen(GenGenInd,10)); % (Pmax-Pmim);
% genRange = [60 60 25 25]'; % HARD-CODED: equal to Pmax-Pmim; where Pmin is 0

% vRange
AllGenInd = mpc.gen(:,1); % 
vRange = mpc.bus(AllGenInd,12) - mpc.bus(AllGenInd,13); %Vmax - Vmin (INCLUDING slack node)

% all variables range
% gvRange = vertcat(genRange,vRange);
% optvar2pu = gvRange;
optvar2pu = genRange./mpc.baseMVA; % Originally for PG only
P_demand_i = sum(mpc.bus(:,PD))/mpc.baseMVA;

% Separating PG and VG variables`
NGen = size(AllGenInd,1);
u_PG = u_NN(1:NGen-1);
u_VG = u_NN(NGen:end);

% Range of the slack bus
SlackgenBusInd = mpc.bus(mpc.bus(:,2) == 3);
[~,SlackgenGenInd] = intersect(mpc.gen(:,1),SlackgenBusInd);
Pslack_max = mpc.gen(SlackgenGenInd,9)/mpc.baseMVA;
Pslack_min = mpc.gen(SlackgenGenInd,10)/mpc.baseMVA;

%% ----- set up iterators -----
% nPoints = size(spaceGrid,1) * size(epsArray,2)

n = 5;
pMin = 5; % $/MW
pMax = 20; % $/MW
% n = 5;
% pMin = 0; % $/MW
% pMax = 5; % $/MW


pVec = linspace(pMin,pMax,n); % COMMENT: LinSpace with evenly distributed prices, for 68 bus system, we may need ccreate it stochastically
% ------ to change only 1 generators prices (eg pslack)
% <<<<<< comment next 2
% [pG1,pG2,pG3] = ndgrid(pVec,pVec,pVec); % meshgrid % 
% spaceGrid = [pG1(:) pG2(:) pG3(:)];
% using lhsdesign
spaceGrid = lhsdesign(length(pVec)^3,3); %%%%%%%%%%% Q: why in power 3 and dimension 3? Should not be dimension n?
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

%% Uniform random samples for gens (alternative to spaceGrid)
NCostSamples = 125;
RandomCostFnc = rand(NCostSamples,NGen);

%% ----- START ITERATION -----
tic
mpcOrigin = mpc;
fprintf('progress:  ')
chr = fprintf('%2.2f %%\n',0);
toc
nk = 1;
% for k = 1:nEps:nPoints
% k = 1;
% for k = 1:NCostSamples   
    
%     disp(['OUTER ITERATION OVER K: ', num2str(k)])
%     obj = spaceGrid(nk,2:end)*(u_PG.*genRange)+spaceGrid(nk,1)*Pslack*mpc.baseMVA;
        obj = RandomCostFnc(nk,2:end)*(u_PG.*genRange) + RandomCostFnc(nk,1)*Pslack*mpc.baseMVA;

    % ----- loss approximation -----
    
    [l,pGeneration] = milp.getLossSensitivity(mpc,mpcOpt,DELTA_P,solverOutput{:}); %% Q: sensitivity for V; how to connect loss with V?
    [l_v,vGeneration] = milp.getLossSensitivity_v(mpc,mpcOpt,DELTA_V,solverOutput{:}); %% Q: sensitivity for V; how to connect loss with V?
    
    for j = 0:nEps-1
%         disp(['INNER ITERATION OVER J: ', num2str(j)])
        i_loss = 0;
        mpc = mpcOrigin;
        err = 1;
        %% LOSS ESTIMATION ITERATION =========================================
        while err >= 0.01  
            % ----- constraints -----
            constr = [];
            constr = [constr;...
                0 <= u_NN <= 1]:'input restriction';%%%%%%%%%%%%% Q: or u_PG?
            constr = [constr;...
                (Pslack_min <= Pslack <= Pslack_max):'slack bus']; %#ok
            % ----- losses -----
            constr = [constr;...
            (Pslack + sum(u_PG.*optvar2pu) - P_demand_i - l(1)...
            + sum(-l(2:NGen).*(u_PG.*optvar2pu - pGeneration(GenGenInd)))... 
            + sum(-l_v(2:NGen+1).*(u_VG.*vRange - vGeneration)) == 0):'lossy power balance'];%#ok
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
            diagnostics = optimize(constr,obj,options);
            milpRes.setpoints  = [value(Pslack)*mpc.baseMVA; value(u_PG).*genRange; value(u_VG).*vRange];
            milpRes.mpc = mpc;

            [err,~,pfSol] = getLosses(milpRes, mpcOpt,NGen);
            fprintf('Loss estimation error: [%d] - %2.3f\n',i_loss,err);
            
            %Added for print by Ilgiz
            milpRes.setpoints

            [l,pGeneration] = iterateLosses(mpc, mpcOpt, pfSol,DELTA_P);
            [l_v,vGeneration] = iterateLosses_VG(mpc, mpcOpt, pfSol,DELTA_V);
            
            mpc = pfSol;
            i_loss = i_loss + 1;
            toc
        end
        % ----- results -----
        setpoints = value(u_PG).*genRange;
        % ----- check solution ----- (will be done in OOPSAT)
%         [milpRes,acopfRes,N1Res] = util.compareResults(mpcOrigin,mpcOpt,obj,setpoints,spaceGrid(nk,:),'silent');
%         [milpRes,acopfRes,N1Res] = util.compareResults(mpcOrigin,mpcOpt,obj,setpoints,RandomCostFnc(nk,:),'silent');
%         pslack(k+j) = value(Pslack);
%         price(k+j) = value(obj);
%         stability(k+j) = milpRes.stability;
%         dampingRatio(k+j) = milpRes.dr; %keep this for backwards compatibility of plots
%         milpRes.time = diagnostics.solvertime;
%         milpRes.diag = diagnostics;
%         gsetpoints(k+j,:) = value(u_PG);%setpoints(:);
%         resultStruct{k+j,1} = milpRes;
%         resultStruct{k+j,2} = acopfRes;
%         resultStruct{k+j,3} = N1Res;
%         
        fprintf('%2.2f %%\n',(k+j)/nPoints*100)
    end
%     nk = nk +1;
% end


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
    l(1) = l(1)/mpc.baseMVA; % HARD-CODED: assuming the slack bus is number 1
end

function [l_v,vGeneration] = iterateLosses_VG(mpc,mpopt,pfResults,DELTA_V)

PG = 2;
PD = 3;
VG = 6;

    pGeneration = pfResults.gen(:, PG);
    vGeneration = pfResults.gen(:, VG);

    genlist = util.GetGenList(mpc); % 0 for slack 1 for rest
    NG = length(genlist);
    l_v = nan(NG,1);

    %% 1 st step -> calculate Pi_DC
    l_v(1) = sum(pGeneration) - sum(mpc.bus(:, PD)); % PD = 3 real power demand

    %% 2 nd step -> ACOPF to calculate P_L^0, \gamma
    for i = 1:NG
        mpcTemp = mpc;
        mpcTemp.gen(:, VG) = vGeneration;
        mpcTemp.gen(i, VG) = vGeneration(i) + DELTA_V;

        pf = runpf(mpcTemp, mpopt);
        loss = sum(pf.gen(:, PG)) - sum(mpc.bus(:, PD));
        l_v(i+1) = (loss - l_v(1))/(DELTA_V*mpc.baseMVA); %divide this one with baseMVA, as delta V is given in p.u.
    end
    
    l_v(1) = l_v(1)/mpc.baseMVA;
end

function [err, varargout] = getLosses(milpRes,mpOpt,NGen)
% get loss comparision of MILP and PF results
% return loss estimation error
% optionally return loss values

    PD = 3; % column numbers, general coded
    PG = 2;
    VG = 6;
    PMAX = 9;
    PMIN = 10;
    VMAX = 12;
    VMIN = 13;
    TOL = 1e-2;
    % References: https://matpower.org/docs/ref/matpower7.0/lib/get_losses.html
    %% Prepare AC PF to run with MILP setpoints
    milpRes.mpc.gen(1:NGen,PMAX) = milpRes.setpoints(1:NGen) + TOL;
    milpRes.mpc.gen(1:NGen,PMIN) = milpRes.setpoints(1:NGen) - TOL;
    milpRes.mpc.gen(1:NGen,PG) = milpRes.setpoints(1:NGen);
    
    milpRes.mpc.gen(1:NGen,VMAX) = milpRes.setpoints(NGen+1:end) + TOL;
    milpRes.mpc.gen(1:NGen,VMIN) = milpRes.setpoints(NGen+1:end) - TOL;
    milpRes.mpc.gen(1:NGen,VG) = milpRes.setpoints(NGen+1:end);
    
%     milpRes.setpoints 
    pfRes = runpf(milpRes.mpc,mpOpt);

    loss.ac_real = sum(real(get_losses(pfRes)));
%     disp(['loss.ac_real: ', num2str(loss.ac_real)])
    
    loss.milp_real = sum(milpRes.setpoints) - sum(milpRes.mpc.bus(:,PD));
%     disp(['loss.milp_real: ', num2str(loss.milp_real)])
    
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
