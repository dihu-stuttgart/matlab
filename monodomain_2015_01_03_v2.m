%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = monodomain_2015_01_03_v2()


%% THIS HAS NEVER BEEN USED FOR THE PHD THESIS AND IS COPIED FROM "COMPUATATIONAL METHODS IN BIOMECHANICS" ASSINGMENT 3 %%% 



  %------------------------------------------------------------------
  % SET PARAMETERS
  % number of grid points in each direction
  n=10;
  % total number of grid points
  num_of_points = n*n;
  % surface area to volume ratio
  A_m = 500.0; % /cm
  % membrane capacitance
  C_m = 1.0; % microF/cm^2
  % effective conductivity
  sigma_eff = 0.01; %3.828; % mS/cm 
  % time step for the PDE
  time_step = 0.1; % ms
  % stop time
  t_end = 10.0; % ms
  % start time of stimulation
  t_start_stimulation = 0.2; % ms
  % stop time of stimulation
  t_end_stimulation = 0.7; % ms
  % grid spacing
  dx = 0.00625; % cm
  % output node number
  out_node = 13;
  
  %------------------------------------------------------------------
  % INITIALISE VARIABLES
  % number of time steps for dynamic PDE solver
  num_of_dt = t_end/time_step;
  % the time step number at which the stimulation starts
  start_stim = t_start_stimulation/time_step;
  % the time step number at which the stimulation stops
  stop_stim = t_end_stimulation/time_step;

  % stimulation points
%  stim_pnts = zeros(n,1);
%   j=2;
%   for i=1:n
%     stim_pnts(i,1) = j;
%     j=j+n;
%   end
  % here, only node 19 is stimulated!!!
  stim_pnts(1,1) = 19;
  % stimulation curtent
  i_Stim = zeros(num_of_points,num_of_dt); 
  
  % membrane voltage as derived from the Hodgkin-Huxley model -- V*
  vS_hh = zeros(num_of_points,1);
  % transmembrane voltage after parabolic PDE evaluation -- V^{k+1}
  V_m = -75*ones(num_of_points,1);
  % transmembrane voltage at the output node for each time step
  V_m_time = zeros(num_of_dt,1);

  % stiffness matrix
  KK = zeros(num_of_points,num_of_points);
  % right hand side vector
  bb  = zeros(num_of_points,1);

  %------------------------------------------------------------------
  % STIMULATION
  % set the stimulation at the specified nodes and times
  for i=start_stim:stop_stim
    for j=1:size(stim_pnts,1)
      node = stim_pnts(j);
      i_Stim(node,i) = 200;%200;
    end
  end
  
  %------------------------------------------------------------------
  % STIFFNESS MATRIX
  % first order dynamic problem -->  KK * V_m = bb
  % introduce factor for compact notation
  fact = time_step/dx^2 * sigma_eff/A_m/C_m;
  for i=n+1:n*n-n
    if i>n
      KK(i,i-n) = -fact;
    end
    KK(i,i-1) = -fact;
    KK(i,i) =  (1+4*fact);
    KK(i,i+1) = -fact;
    if i<n*n-n+1
      KK(i,i+n) = -fact;
    end
  end
  
  % BOUNDARY CONDITIONS
  % homogeneous Neumann B.C. - no flow
  % bottom side
  for i=1:n
    KK(i,:) = 0;
    KK(i,i) = 1;
    KK(i,i+n) = -1;
  end
  % top side
  for i=n*n-n+1:n*n
    KK(i,:) = 0;
    KK(i,i) = 1;
    KK(i,i-n) = -1;
  end
  % right side
  for i=2*n:n:n*n-n
    KK(i,:) = 0;
    KK(i,i) = 1;
    KK(i,i-1) = -1;
  end
  % left side
  for i=n+1:n:n*n-2*n+1
    KK(i,:) = 0;
    KK(i,i) = 1;
    KK(i,i+1) = -1;
  end
  
  
  %---------------------------------------------------------------------
  % CELLULAR MODEL
  
  % variables in cellular model:
  % VOI       -- time
  % CONSTANTS -- constants
  % ALGEBRAIC -- algebraic variables, e.g. i_Na, i_K
  %              (no time derivative exists for algebraic variables)
  % STATES    -- state variables, e.g. V
  % RATES     -- time derivative of state variables, e.g. d(V)/dt
  
  % number of algebraic variables in the system
  global algebraicVariableCount;  
  algebraicVariableCount = getAlgebraicVariableCount();

  % Set numerical accuracy options for ODE solver
  options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 0.1);
  
  % Initialise constants and state variables for cellular model
  [CONSTANTS] = initConsts();
  [INIT_STATES] = initStates(); 
  
  % store all state variables for all points
  ALL_STATES = zeros(num_of_points,4); % 4 is the number of STATE variables
  % initialise ALL_STATES
  for i = 1:num_of_points
    ALL_STATES(i,:) = INIT_STATES;
  end
  
  %---------------------------------------------------------------------
  % SOLUTION PROCESS
  %---------------------------------------------------------------------
  
  % loop over the PDE time steps
  for time = 1:num_of_dt
    
    fprintf('step #   : %d \t', time);
    fprintf('time [ms]: %f \n', (time-1)*time_step);
    % Set timespan to integrate over 
    % NOTE: this does not affect the time step for the ODE/DAE integration
    % but only specifies the times for the output. 
    % [first last] returns values at all integration time steps!!!
    tspan = [(time-1)*time_step, (time-0.5)*time_step, time*time_step];

    %---------------------------------------------------------------------
    % CELLULAR MODEL
    %---------------------------------------------------------------------
    
    % Integrate the cellular model at each discretisation point
    for point = 1:num_of_points
      
      % load the last state as initial state for the next integration
      LAST_STATES = ALL_STATES(point,:);
      
      %---------------------------------------------------------------------
      % INTEGRATE CELLULAR MODEL WITH ODE/DAE SOLVER
      % NOTE: no timestep for the integration of the cellular model has to
      % be provided. MATLAB chooses an appropriate time step
      
      [VOI, STATES] = ode15s(@(VOI, STATES)computeRates(VOI, STATES, CONSTANTS, i_Stim(point,time)), tspan, LAST_STATES, options);
      
      % Compute algebraic variables
%      [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS, I_HH(point,time));
%      ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI, I_HH(point,time));
      
      % update ALL_STATES
      ALL_STATES(point,:) = STATES(end,:);
      % update transmembrane voltage of cellular model
      vS_hh(point,1) = STATES(end,1);
      
    end % points

    % update the transmembrane voltage V_m
    V_m = vS_hh(:,1);
    
    %---------------------------------------------------------------------
    % PARABOLIC EQUATION
    %---------------------------------------------------------------------

    % set the right hand side vector
    bb(:,1) = V_m(:,1);
    
    %---------------------------------------------------------------------
    % BOUNDARY CONDITIONS
    % homogeneous Neumann B.C., no flow 
    % bottom side
    for i=1:n
      bb(i,1) = 0;
    end
    % top side
    for i=n*n-n+1:n*n
      bb(i,1) = 0;
    end
    % right side
    for i=2*n:n:n*n-n
      bb(i,1) = 0;
    end
    % left side
    for i=n+1:n:n*n-2*n+1
      bb(i,1) = 0;
    end
    
    %---------------------------------------------------------------------
    % SOLVE THE PARABOLIC PDE
    V_m = KK\bb;
    
    %---------------------------------------------------------------------
    % update the cell models transmembrane voltage
    ALL_STATES(:,1)  = V_m;
    % store the transmembrane voltage at the output node
    V_m_time(time,1) = V_m(out_node,1);
    
    % TODO tomo
%    x_a = linspace(0,0.025,n);
    y_a = (0:1/(n-1):1);
    x_a = y_a;
    
    % surface plots for time = 1,2,3,.. ms
    if (mod(time,1)==0)
      % change format of solution from vector to matrix
      u = zeros(n,n);
      for i = 1:n
        u(:,i) = V_m((i-1)*n+1:i*n,1);
      end;
      % plotting
      figure(time);
      surf(meshgrid(x_a)',meshgrid(y_a),u);
      axis([0 1 0 1 -80 40 0 1]);
    end

  end %time
  
  tt = linspace(time_step,t_end,num_of_dt);
  figure(99);
  plot(tt, V_m_time);
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                STRICTLY NO CHANGES BELOW THIS POINT                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% THE CELLULAR MODDEL OF HODGKIN & HUXLEY, 1952

function [algebraicVariableCount] = getAlgebraicVariableCount() 
    % Used later when setting a global variable with the number of algebraic variables.
    % Note: This is not the "main method".  
    algebraicVariableCount =10;
end
% There are a total of 4 entries in each of the rate and state variable arrays.
% There are a total of 8 entries in the constant variable array.

function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_VOI = ''; LEGEND_CONSTANTS = '';
    LEGEND_VOI = strpad('time in component environment (millisecond)');
    LEGEND_STATES(:,1) = strpad('V in component membrane (millivolt)');
    LEGEND_CONSTANTS(:,1) = strpad('E_R in component membrane (millivolt)');
    LEGEND_CONSTANTS(:,2) = strpad('Cm in component membrane (microF_per_cm2)');
    LEGEND_ALGEBRAIC(:,5) = strpad('i_Na in component sodium_channel (microA_per_cm2)');
    LEGEND_ALGEBRAIC(:,9) = strpad('i_K in component potassium_channel (microA_per_cm2)');
    LEGEND_ALGEBRAIC(:,10) = strpad('i_L in component leakage_current (microA_per_cm2)');
    LEGEND_ALGEBRAIC(:,1) = strpad('i_Stim in component membrane (microA_per_cm2)');
    LEGEND_CONSTANTS(:,3) = strpad('g_Na in component sodium_channel (milliS_per_cm2)');
    LEGEND_CONSTANTS(:,6) = strpad('E_Na in component sodium_channel (millivolt)');
    LEGEND_STATES(:,2) = strpad('m in component sodium_channel_m_gate (dimensionless)');
    LEGEND_STATES(:,3) = strpad('h in component sodium_channel_h_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,2) = strpad('alpha_m in component sodium_channel_m_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,6) = strpad('beta_m in component sodium_channel_m_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,3) = strpad('alpha_h in component sodium_channel_h_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,7) = strpad('beta_h in component sodium_channel_h_gate (per_millisecond)');
    LEGEND_CONSTANTS(:,4) = strpad('g_K in component potassium_channel (milliS_per_cm2)');
    LEGEND_CONSTANTS(:,7) = strpad('E_K in component potassium_channel (millivolt)');
    LEGEND_STATES(:,4) = strpad('n in component potassium_channel_n_gate (dimensionless)');
    LEGEND_ALGEBRAIC(:,4) = strpad('alpha_n in component potassium_channel_n_gate (per_millisecond)');
    LEGEND_ALGEBRAIC(:,8) = strpad('beta_n in component potassium_channel_n_gate (per_millisecond)');
    LEGEND_CONSTANTS(:,5) = strpad('g_L in component leakage_current (milliS_per_cm2)');
    LEGEND_CONSTANTS(:,8) = strpad('E_L in component leakage_current (millivolt)');
    LEGEND_RATES(:,1) = strpad('d/dt V in component membrane (millivolt)');
    LEGEND_RATES(:,2) = strpad('d/dt m in component sodium_channel_m_gate (dimensionless)');
    LEGEND_RATES(:,3) = strpad('d/dt h in component sodium_channel_h_gate (dimensionless)');
    LEGEND_RATES(:,4) = strpad('d/dt n in component potassium_channel_n_gate (dimensionless)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end

function [CONSTANTS] = initConsts()
    VOI = 0; CONSTANTS = [];
    CONSTANTS(:,1) = -75;
    CONSTANTS(:,2) = 1;
    CONSTANTS(:,3) = 120;
    CONSTANTS(:,4) = 36;
    CONSTANTS(:,5) = 0.3;
    CONSTANTS(:,6) = CONSTANTS(:,1)+115.000;
    CONSTANTS(:,7) = CONSTANTS(:,1) - 12.0000;
    CONSTANTS(:,8) = CONSTANTS(:,1)+10.6130;
end

function [STATES] = initStates()
    STATES(:,1) = -75;
    STATES(:,2) = 0.05;
    STATES(:,3) = 0.6;
    STATES(:,4) = 0.325;
    if (isempty(STATES)), warning('Initial values for states not set'); end
end

function [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS, i_Stim)
    global algebraicVariableCount;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        ALGEBRAIC = zeros(1, algebraicVariableCount);
    else
        statesRowCount = statesSize(1);
        ALGEBRAIC = zeros(statesRowCount, algebraicVariableCount);
        RATES = zeros(statesRowCount, statesColumnCount);
    end
    ALGEBRAIC(:,2) = (  - 0.100000.*(STATES(:,1)+50.0000))./((exp(( - (STATES(:,1)+50.0000)./10.0000))) - 1.00000);
    ALGEBRAIC(:,6) =  4.00000.*(exp(( - (STATES(:,1)+75.0000)./18.0000)));
    RATES(:,2) =  ALGEBRAIC(:,2).*(1.00000 - STATES(:,2)) -  ALGEBRAIC(:,6).*STATES(:,2);
    ALGEBRAIC(:,3) =  0.0700000.*(exp(( - (STATES(:,1)+75.0000)./20.0000)));
    ALGEBRAIC(:,7) = 1.00000./((exp(( - (STATES(:,1)+45.0000)./10.0000)))+1.00000);
    RATES(:,3) =  ALGEBRAIC(:,3).*(1.00000 - STATES(:,3)) -  ALGEBRAIC(:,7).*STATES(:,3);
    ALGEBRAIC(:,4) = (  - 0.0100000.*(STATES(:,1)+65.0000))./((exp(( - (STATES(:,1)+65.0000)./10.0000))) - 1.00000);
    ALGEBRAIC(:,8) =  0.125000.*(exp(((STATES(:,1)+75.0000)./80.0000)));
    RATES(:,4) =  ALGEBRAIC(:,4).*(1.00000 - STATES(:,4)) -  ALGEBRAIC(:,8).*STATES(:,4);
    ALGEBRAIC(:,5) =  CONSTANTS(:,3).*(STATES(:,2) .^ 3.00000).*STATES(:,3).*(STATES(:,1) - CONSTANTS(:,6));
    ALGEBRAIC(:,9) =  CONSTANTS(:,4).*(STATES(:,4) .^ 4.00000).*(STATES(:,1) - CONSTANTS(:,7));
    ALGEBRAIC(:,10) =  CONSTANTS(:,5).*(STATES(:,1) - CONSTANTS(:,8));
% tomo
%    ALGEBRAIC(:,1) = piecewise({VOI>=10.0000&VOI<=10.5000, 20.0000 }, 0.00000);
    ALGEBRAIC(:,1) = i_Stim;
% tomo end
    RATES(:,1) =  - ( - ALGEBRAIC(:,1)+ALGEBRAIC(:,5)+ALGEBRAIC(:,9)+ALGEBRAIC(:,10))./CONSTANTS(:,2);
   RATES = RATES';
end

% Calculate algebraic variables
function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI, i_Stim)
    ALGEBRAIC(:,2) = (  - 0.100000.*(STATES(:,1)+50.0000))./((exp(( - (STATES(:,1)+50.0000)./10.0000))) - 1.00000);
    ALGEBRAIC(:,6) =  4.00000.*(exp(( - (STATES(:,1)+75.0000)./18.0000)));
    ALGEBRAIC(:,3) =  0.0700000.*(exp(( - (STATES(:,1)+75.0000)./20.0000)));
    ALGEBRAIC(:,7) = 1.00000./((exp(( - (STATES(:,1)+45.0000)./10.0000)))+1.00000);
    ALGEBRAIC(:,4) = (  - 0.0100000.*(STATES(:,1)+65.0000))./((exp(( - (STATES(:,1)+65.0000)./10.0000))) - 1.00000);
    ALGEBRAIC(:,8) =  0.125000.*(exp(((STATES(:,1)+75.0000)./80.0000)));
    ALGEBRAIC(:,5) =  CONSTANTS(:,3).*(STATES(:,2) .^ 3.00000).*STATES(:,3).*(STATES(:,1) - CONSTANTS(:,6));
    ALGEBRAIC(:,9) =  CONSTANTS(:,4).*(STATES(:,4) .^ 4.00000).*(STATES(:,1) - CONSTANTS(:,7));
    ALGEBRAIC(:,10) =  CONSTANTS(:,5).*(STATES(:,1) - CONSTANTS(:,8));
% tomo
%    ALGEBRAIC(:,1) = piecewise({VOI>=10.0000&VOI<=10.5000, 20.0000 }, 0.00000);
    ALGEBRAIC(:,1) = i_Stim;
% tomo end
end

% Compute result of a piecewise function
function x = piecewise(cases, default)
    set = [0];
    for i = 1:2:length(cases)
        if (length(cases{i+1}) == 1)
            x(cases{i} & ~set,:) = cases{i+1};
        else
            x(cases{i} & ~set,:) = cases{i+1}(cases{i} & ~set);
        end
        set = set | cases{i};
        if(set), break, end
    end
    if (length(default) == 1)
        x(~set,:) = default;
    else
        x(~set,:) = default(~set);
    end
end

% Pad out or shorten strings to a set length
function strout = strpad(strin)
    req_length = 160;
    insize = size(strin,2);
    if insize > req_length
        strout = strin(1:req_length);
    else
        strout = [strin, blanks(req_length - insize)];
    end
end

