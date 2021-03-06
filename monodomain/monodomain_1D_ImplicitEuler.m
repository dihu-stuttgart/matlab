%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = monodomain_1D_ImplicitEuler()


%% THIS HAS NEVER BEEN USED FOR THE PHD THESIS AND IS COPIED FROM "COMPUATATIONAL METHODS IN BIOMECHANICS" ASSINGMENT 3 %%% 
  clc;
  close all; 
  clear all;
  
  %------------------------------------------------------------------
  % SET PARAMETERS
  % number of grid points in each direction
  n=24;%default 10
  % total number of grid points
  num_of_points = n;
  
  % Length
  L=1; %cm
  % grid spacing
  dx = L/n; % cm
  
  % surface area to volume ratio
  A_m = 500.0; % /cm
  % membrane capacitance, fast-twitch
  C_m = 1.0; % microF/cm^2
  % effective conductivity
  sigma_eff = 3.828; %3.828; % mS/cm 
  % Diffusion coeff.
  D=sigma_eff/A_m/C_m;
  
  % time step for the PDE
  time_step_pde = 0.05;
% time step for the ODE
 % time_step_ode = 0.05;

  % stop time
  t_end = 5.0; % ms
  
  %frequency of stimulation
  %f=60;
  
  % start time of stimulation
  t_start_stimulation = 0.0;
  % stop time of stimulation
  t_end_stimulation = 0.1;
  
  %fast twitch
  I_stim=2000;
  
  %output time
  t_out=3;%ms  
  
  %------------------------------------------------------------------
  % INITIALISE VARIABLES
  % number of time steps for dynamic PDE solver
  num_of_dt = t_end/time_step_pde;
  % the time step number at which the stimulation starts
  if(t_start_stimulation ~= 0.0) 
    start_stim = t_start_stimulation/time_step_pde;
  else
    start_stim=1;
  end
  % the time step number at which the stimulation stops
  stop_stim = t_end_stimulation/time_step_pde;

  % stimulation points
  %stim_pnts = zeros(n);
  % here, only one node in the middle of the domain is stimulated!!!
  stim_pnts(1) = ceil(n/2.0);
  % stimulation curtent
  i_Stim = zeros(num_of_points,num_of_dt); 
  
  % membrane voltage as derived from the Hodgkin-Huxley model -- V*
  vS_hh = zeros(num_of_points,1);
  % transmembrane voltage after parabolic PDE evaluation -- V^{k+1}
  V_m = -75*ones(num_of_points,1);
  % transmembrane voltage at the output node for each time step
  V_m_time = zeros(num_of_dt,num_of_points);

  %------------------------------------------------------------------
  % STIMULATION
  % set the stimulation at the specified nodes and times
  for i=start_stim:stop_stim
    for j=1:length(stim_pnts)
      node = stim_pnts(j);
      i_Stim(node,i) = I_stim;%1 not enough stimulation
    end
  end
  
  %-----------------------------------------------------------------
  %method of time discretization
  method='ImplicitEuler';
  
  %------------------------------------------------------------------
  % STIFFNESS MATRIX
  % first order dynamic problem -->  KK * V_m = bb
  KK= StiffnessMatrix(method,num_of_points,time_step_pde,dx,D);
  
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
  algebraicVariableCount = getAlgebraicVariableCount('HODGKIN_HUXLEY');

  % Set numerical accuracy options for ODE solver
  options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', time_step_pde);
  
  % Initialise constants and state variables for cellular model
  [CONSTANTS] = initConsts('HODGKIN_HUXLEY');
  [INIT_STATES] = initStates('HODGKIN_HUXLEY'); 
  
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
    
    %fprintf('step #   : %d \t', time);
    %fprintf('time [ms]: %f \n', (time-1)*time_step);
    % Set timespan to integrate over 
    % NOTE: this does not affect the time step for the ODE/DAE integration
    % but only specifies the times for the output. 
    % [first last] returns values at all integration time steps!!!
    tspan = [(time-1)*time_step_pde, (time-0.5)*time_step_pde, time*time_step_pde];

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
      
      [VOI, STATES] = ode15s(@(VOI, STATES)computeRates_HODGKIN_HUXLEY(VOI, STATES, CONSTANTS, i_Stim(point,time)), tspan, LAST_STATES, options);
%       if (point==1)
%       STATES
%       end
      
      % Compute algebraic variables
%      [RATES, ALGEBRAIC] = computeRates_HODGKIN_HUXLEY(VOI, STATES, CONSTANTS, i_stim(point,time));
%      ALGEBRAIC = computeAlgebraic_HODGKIN_HUXLEY(ALGEBRAIC, CONSTANTS, STATES, VOI, I_HH(point,time));
      
      % update ALL_STATES
      ALL_STATES(point,:) = STATES(end,:);
      % update transmembrane voltage of cellular model
      vS_hh(point,1) = STATES(end,1);
      
    end % points

    % update the transmembrane voltage V_m
    V_m = vS_hh;
    
    %---------------------------------------------------------------------
    % PARABOLIC EQUATION
    %---------------------------------------------------------------------
     bb=setRHS(method,V_m);
        
    %---------------------------------------------------------------------
    % SOLVE THE PARABOLIC PDE
    V_m = KK\bb;
    
    %---------------------------------------------------------------------
    % update the cell models transmembrane voltage
    ALL_STATES(:,1)  = V_m;
    % store the transmembrane voltage at the output node
    V_m_time(time,:) = V_m;

  end %time
  
  OutToFile(V_m_time,t_out,time_step_pde,method); 
  
  tt = linspace(0,t_end,num_of_dt);
  x_a = linspace(0,1,num_of_points);
  surf(x_a,tt,V_m_time);
  
end

function []=OutToFile(V_m_time,t_out,time_step,method)
time=t_out/time_step;

outfile=fopen(strcat('out','Vm_',num2str(time_step),'_',method,'.txt'),'w');
fprintf(outfile,'t_end=%f\ntime_step=%f\n',t_out,time_step);
fprintf(outfile,'%6.6f\n',V_m_time(time,:));
end
