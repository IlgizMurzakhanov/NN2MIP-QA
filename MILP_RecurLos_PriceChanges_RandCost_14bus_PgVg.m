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
epsArray = 10;
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
% mpc = util.case14_wind;     % base mpc case
mpc = util.case14_wind_corr;     % base mpc case

%% modify mpc (Only for 14 bus system)
% % ====== relax q limits
% QMAX = 4;
% QMIN = 5;
% mpc.gen(:,QMAX) = mpc.gen(:,QMAX)+0.25*(mpc.gen(:,QMAX)-mpc.gen(:,QMIN));
% mpc.gen(:,QMIN) = mpc.gen(:,QMIN)-0.25*(mpc.gen(:,QMAX)-mpc.gen(:,QMIN));
% % ===== Relax S limits
% sTol = 1.3;
% mpc.branch(:,6) = mpc.branch(:,6)*sTol; 
% mpc.branch(:,7) = mpc.branch(:,7)*sTol;
% mpc.branch(:,8) = mpc.branch(:,7)*sTol;

%% mpopt
mpcOpt = mpoption;
mpcOpt.verbose = 0;
mpcOpt.out.all = 0;

qOpt = 1;
mpcOpt.pf.enforce_q_lims = qOpt;
solverOutput = {
    'verbose',0};

%%%%%%%%%%%%%%%%%%%%%% SET THE PATH: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path = [pwd,filesep,'neural_network',filesep,'qLimitsNotEnforcedQS3_2relaxed',filesep];
% path = 'C:\Users\??????\Documents\?????????_??????\DTU\Research\Collaboration with Andreas\neural_network\14_bus_system_VG\';
% path = 'C:\Users\??????\Documents\?????????_??????\DTU\Research\Collaboration with Andreas\neural_network\14_bus_system_PG_VG\dataset_181120\';
path = 'C:\Users\??????\Documents\?????????_??????\DTU\Research\Collaboration with Andreas\neural_network\Dataset_14bus_PG_VG\';

HIDDEN_LAYER = 50;          % MILP formulation
RELU_LAYERS = 3;            % MILP formulation
DELTA_P = 10^-1;               % losses % HARD-CODED: see paper, given in MW
DELTA_V = 10^-1;             % given in p.u.

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

%% mpc case data preprocessing
% genRange
SlackBusInd = mpc.bus(mpc.bus(:,2) == 3); % slack bus index
[~,SlackGenInd] = intersect(mpc.gen(:,1),SlackBusInd); % .. in gen table ..

GenBusInd = mpc.bus(mpc.bus(:,2) == 2); % indices of gen-s in bus table (EXCLUDING slack node)
[~,GenGenInd] = intersect(mpc.gen(:,1),GenBusInd); % .. in gen table ..

genRange = (mpc.gen(GenGenInd,9) - mpc.gen(GenGenInd,10)); % (Pmax-Pmim);
optvar2pu = genRange./mpc.baseMVA; % Originally for PG only
P_demand_i = sum(mpc.bus(:,PD))/mpc.baseMVA;

% voltage range
AllGenBusInd = mpc.gen(:,1); % 
Vmax = mpc.bus(AllGenBusInd,12);
Vmin = mpc.bus(AllGenBusInd,13);
vRange = Vmax - Vmin; %Vmax - Vmin (INCLUDING slack node)

% Number of generators
NGenIncSlack = size(AllGenBusInd,1); % including the slack node
NGenExcSlack = size(GenBusInd,1); % excluding ..

% Separating PG and VG variables (is not used further)
% u_PG = u_NN(1:NGenExcSlack); % All generators except slack are going first
% u_VG = u_NN(NGenIncSlack:end); % The rest are voltages, on generator buses, including slack

% Range of the slack bus
SlackgenBusInd = mpc.bus(mpc.bus(:,2) == 3);
[~,SlackgenGenInd] = intersect(mpc.gen(:,1),SlackgenBusInd);
Pslack_max = mpc.gen(SlackgenGenInd,9)/mpc.baseMVA;
Pslack_min = mpc.gen(SlackgenGenInd,10)/mpc.baseMVA;

%% ----- set up iterators -----
% % nPoints = size(spaceGrid,1) * size(epsArray,2)
% 
% n = 1;
% pMin = 5; % $/MW
% pMax = 5; % $/MW
% % n = 5;
% % pMin = 0; % $/MW
% % pMax = 5; % $/MW
% 
% 
% pVec = linspace(pMin,pMax,n); % COMMENT: LinSpace with evenly distributed prices, for 68 bus system, we may need ccreate it stochastically
% % ------ to change only 1 generators prices (eg pslack)
% % <<<<<< comment next 2
% % [pG1,pG2,pG3] = ndgrid(pVec,pVec,pVec); % meshgrid % 
% % spaceGrid = [pG1(:) pG2(:) pG3(:)];
% 
% % using lhsdesign
% spaceGrid = lhsdesign(length(pVec)^3,3); %%%%%%%%%%% Q: why in power 3 and dimension 3? Should not be dimension n?
% spaceGrid = (spaceGrid.*(pMax-pMin))+pMin;
% 
% % <<<<<< insert next
% spaceGrid = [pVec', repmat(mpc.gencost(2,6),[n,1]),repmat(mpc.gencost(3,6),[n,1])];
% 
% nEps = size(epsArray,2);
% nCostGrid = size(spaceGrid,1);
% nPoints = nCostGrid * nEps;
% pslack= nan(nPoints, 1);
% price = nan(nPoints, 1); % price (for coloring)
% stability = nan(nPoints, 1); % see if the solution did not violate constraints
% dampingRatio = stability;
% gsetpoints = nan(nPoints, 4);
% vsetpoints = nan(nPoints, 5);
% resultStruct = cell(nPoints,3); %milp,acopf, n-1 acopf
% 
% % FIX generator 6 and 8 cost to 0
% spaceGrid = [spaceGrid zeros(nCostGrid,2) ];

%% Uniform random samples for gens (alternative to spaceGrid)
nEps = size(epsArray,2);
NCostSamples = 5;
RandomCostFnc = rand(NCostSamples,NGenIncSlack);
nPoints = NCostSamples * nEps;

% Repeated from SpaceGrid
pslack= nan(nPoints, 1);
price = nan(nPoints, 1); % price (for coloring)
stability = nan(nPoints, 1); % see if the solution did not violate constraints
dampingRatio = stability;
gsetpoints = nan(nPoints, 4);
vsetpoints = nan(nPoints, 5);
resultStruct = cell(nPoints,3); %milp,acopf, n-1 acopf

%% For analysis of iterative process
NStopIter = 20; % number of stopping iterations in err loop
Output_cell = cell(9,NStopIter);

%% ----- START ITERATION -----
tic
mpcOrigin = mpc;
fprintf('progress:  ')
chr = fprintf('%2.2f %%\n',0);
toc
nk = 1;
for k = 1:nEps:nPoints
    
%     disp(['ITERATION OVER k (nPoints) : ', num2str(k)])

%     %% SpaceGrid
%     obj = spaceGrid(nk,2:end)*((u_NN(1:NGenExcSlack)).*genRange)+spaceGrid(nk,1)*Pslack*mpc.baseMVA;
    
    %% RandomCost
    obj = RandomCostFnc(nk,2:end)*((u_NN(1:NGenExcSlack)).*genRange) + RandomCostFnc(nk,1)*Pslack*mpc.baseMVA;

    %% ----- loss approximation -----
%     disp('Initial loss and pGeneration')
    [l,pGeneration] = milp.getLossSensitivity_p(mpc,mpcOpt,DELTA_P,NGenIncSlack,solverOutput{:}); %% Q: sensitivity for V; how to connect loss with V?
    
%     disp('Initial loss and vGeneration')
    [l_v,vGeneration] = milp.getLossSensitivity_v(mpc,mpcOpt,DELTA_V,SlackBusInd,SlackGenInd,GenBusInd,GenGenInd,NGenIncSlack,solverOutput{:}); %% Q: sensitivity for V; how to connect loss with V?
%     vGeneration_fixed = vGeneration;
%     vGeneration_fixed = ones(NGenIncSlack,1);
%     vGeneration_fixed = Vmax; 
%     vGeneration_fixed = Vmin; 
    vGeneration_fixed = [1.06;1.045;1.01;1.02;1.01]; 
    
    for j = 0:nEps-1
%         disp(['ITERATION OVER j (nEps) : ', num2str(j)])
        i_loss = 0;
        mpc = mpcOrigin;
        err = 1;
        
        %% LOSS ESTIMATION ITERATION =========================================
        while err >= 0.01  
            disp(['Iteration ', num2str(i_loss), ' in the err loop'])
            % ----- constraints -----
            constr = [];
            constr = [constr;...
                0 <= u_NN <= 1]:'input restriction';%%%%%%%%%%%%% Q: or u_PG?
            constr = [constr;...
                (Pslack_min <= Pslack <= Pslack_max):'slack bus']; %#ok
            % ----- losses -----
            constr = [constr;...
            (Pslack + sum((u_NN(1:NGenExcSlack)).*optvar2pu) - P_demand_i - l(1)...
            + sum(-l(2:NGenIncSlack).*((u_NN(1:NGenExcSlack)).*optvar2pu - pGeneration(GenGenInd)))... 
            + sum(-l_v(2:NGenIncSlack+1).*(Vmin + (u_NN(NGenIncSlack:end)).*vRange - vGeneration)) == 0):'lossy power balance'];%#ok
            % ------- voltage fixation ----------
            constr = [constr;...
                (Vmin + (u_NN(NGenIncSlack:end)).*vRange == vGeneration_fixed):'voltage fixation'];%#
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
            diagnostics.info
            milpRes.setpoints  = [value(Pslack)*mpc.baseMVA; value(u_NN(1:NGenExcSlack)).*genRange; Vmin + value(u_NN(NGenIncSlack:end)).*vRange];
%             disp('milpRes.setpoints')
%             milpRes.setpoints
            
            milpRes.mpc = mpc;

            
            [err,PFstatus,~,pfSol] = getLosses(milpRes, mpcOpt,NGenIncSlack);
            fprintf('Loss estimation error (abs(loss.ac_real - loss.milp_real)/loss.ac_real) : [%d] - %2.3f\n',i_loss,err);
            
            [l,pGeneration] = iterateLosses(mpc, mpcOpt, pfSol,DELTA_P,NGenIncSlack);
            [l_v,vGeneration] = iterateLosses_VG(mpc, mpcOpt, pfSol,DELTA_V,SlackBusInd,SlackGenInd,GenBusInd,GenGenInd,NGenIncSlack);
            
            % Filling the cell
            Output_cell{1,i_loss+1} = l;
            Output_cell{2,i_loss+1} = l_v;
            Output_cell{3,i_loss+1} = pGeneration;
            Output_cell{4,i_loss+1} = vGeneration;
            Output_cell{6,i_loss+1} = diagnostics.info;
            Output_cell{5,i_loss+1} = milpRes.setpoints;
            Output_cell{7,i_loss+1} = PFstatus;
            Output_cell{8,i_loss+1} = err;
            Output_cell{9,i_loss+1} = value(obj);
            
            mpc = pfSol;
            i_loss = i_loss + 1;
            if i_loss > NStopIter % it should not take that long
                disp('terminated, as err loop exceeded set number of iterations')
                break
            end
            toc
        end
        % ----- results -----
        setpoints = [value(u_NN(1:NGenExcSlack)).*genRange; Vmin + value(u_NN(NGenIncSlack:end)).*vRange];
        % ----- check solution ----- (will be done in OOPSAT)
%         [milpRes,acopfRes,N1Res] = util.compareResults(mpcOrigin,mpcOpt,obj,setpoints,RandomCostFnc(nk,:),'silent');
        [milpRes,acopfRes,N1Res] = util.compareResults(mpcOrigin,mpcOpt,obj,setpoints,RandomCostFnc(nk,:),'loud');
        pslack(k+j) = value(Pslack);
        price(k+j) = value(obj);
        stability(k+j) = milpRes.stability;
        dampingRatio(k+j) = milpRes.dr; %keep this for backwards compatibility of plots
        milpRes.time = diagnostics.solvertime;
        milpRes.diag = diagnostics;
        gsetpoints(k+j,:) = u_NN(1:NGenExcSlack);%setpoints(:);
        vsetpoints(k+j,:) = u_NN(NGenIncSlack:end);%setpoints(:);
        resultStruct{k+j,1} = milpRes;
        resultStruct{k+j,2} = acopfRes;
        resultStruct{k+j,3} = N1Res;
        
        fprintf('%2.2f %%\n',(k+j)/nPoints*100)
    end
    nk = nk +1;
end

%% Saving file
toc
warning('on')
fileName = strcat(fileNamePre,'.mat');
save(fileName,'setpoints','pslack','price','stability','gsetpoints','vsetpoints','dampingRatio','resultStruct','epsArray')
%%%% saved variables
% pslack        - value of the slack generator real power output
% spaceGrid     - the grid of cost functions (each row defines a cost function)
% price         - the price of the MILP solution
% stability     - 1 if MILP satisfies all criteria
% gsetpoints    - generator setpoints returned by MILP
% vsetpoints    - voltage setpoints returned by MILP
% damping ratio - damping ratio

diary off
if log
    fprintf('Open <a href="matlab: opentoline(''%s'',1)">log file</a>!\n',logFile) %#ok
end

%% ===== END =====








%% ----- functions -----
function [l,pGeneration] = iterateLosses(mpc,mpopt,pfResults,DELTA_P,NGenIncSlack)
PG = 2;
PD = 3;
    pGeneration = pfResults.gen(:, PG);
    l = nan(NGenIncSlack,1);

    %% 1 st step -> calculate Pi_DC
    l(1) = sum(pGeneration) - sum(mpc.bus(:, PD)); % PD = 3 real power demand

    %% 2 nd step -> ACOPF to calculate P_L^0, \gamma
    for i = 2:NGenIncSlack
        mpcTemp = mpc;
        mpcTemp.gen(:, PG) = pGeneration;
        mpcTemp.gen(i, PG) = pGeneration(i)+ DELTA_P;

        pf = runpf(mpcTemp, mpopt);
        loss = sum(pf.gen(:, PG)) - sum(mpc.bus(:, PD));
        l(i) = (loss - l(1))/DELTA_P; % don't divide this one with baseMVA, as delta P is given in MW
    end
    pGeneration = pGeneration./mpc.baseMVA;
    l(1) = l(1)/mpc.baseMVA; % HARD-CODED: assuming the slack bus is number 1
end

function [l_v,vGeneration] = iterateLosses_VG(mpc,mpopt,pfResults,DELTA_V,SlackBusInd,SlackGenInd,GenBusInd,GenGenInd,NGenIncSlack)

    PG = 2;
    PD = 3;
    VG = 6;
    VB = 8;

    pGeneration = pfResults.gen(:, PG);
   
    % For the case when slack bus is not #1
    vGeneration(SlackGenInd) = pfResults.bus(SlackBusInd, VB);
    vGeneration(GenGenInd) = pfResults.bus(GenBusInd, VB);
    vGeneration = vGeneration';
    
    l_v = nan(NGenIncSlack+1,1);

    %% 1 st step -> calculate Pi_DC
    l_v(1) = sum(pGeneration) - sum(mpc.bus(:, PD)); 

    %% 2 nd step -> ACOPF to calculate P_L^0, \gamma
    for i = 1:NGenIncSlack
        mpcTemp = mpc;
        mpcTemp.gen(i, VG) = vGeneration(i) + DELTA_V;

        pf = runpf(mpcTemp, mpopt);
        loss = sum(pf.gen(:, PG)) - sum(mpc.bus(:, PD));
        l_v(i+1) = (loss - l_v(1))/(DELTA_V*mpc.baseMVA); %divide this one with baseMVA, as delta V is given in p.u.
    end
    
    l_v(1) = l_v(1)/mpc.baseMVA; % losses, in p.u.
    
end

function [err, PFstatus, varargout] = getLosses(milpRes,mpOpt,NGenIncSlack)
% get loss comparision of MILP and PF results
% return loss estimation error
% optionally return loss values

    
    PG = 2;
    PD = 3; 
    VG = 6;
    VB = 8;
    
    PMAX = 9;
    PMIN = 10;
    VMAX = 12;
    VMIN = 13;
    
    TOL = 1e-2;
    TOL_V = 1e-3;
    % References: https://matpower.org/docs/ref/matpower7.0/lib/get_losses.html
    
    %% Prepare AC PF to run with MILP setpoints
    milpRes.mpc.gen(1:NGenIncSlack,PMAX) = milpRes.setpoints(1:NGenIncSlack) + TOL;
    milpRes.mpc.gen(1:NGenIncSlack,PMIN) = milpRes.setpoints(1:NGenIncSlack) - TOL;
    milpRes.mpc.gen(1:NGenIncSlack,PG) = milpRes.setpoints(1:NGenIncSlack);
    
    milpRes.mpc.gen(1:NGenIncSlack,VMAX) = milpRes.setpoints(NGenIncSlack+1:end) + TOL_V;
    milpRes.mpc.gen(1:NGenIncSlack,VMIN) = milpRes.setpoints(NGenIncSlack+1:end) - TOL_V;
    milpRes.mpc.gen(1:NGenIncSlack,VG) = milpRes.setpoints(NGenIncSlack+1:end);
    
    pfRes = runpf(milpRes.mpc,mpOpt);
    disp('Power flow execution: success - 1, failure - 0')
    PFstatus = pfRes.success

    loss.ac_real = sum(real(get_losses(pfRes)));
    loss.milp_real = sum(milpRes.setpoints(1:NGenIncSlack)) - sum(milpRes.mpc.bus(:,PD));
    
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
