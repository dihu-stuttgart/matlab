%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = monodomain_1D_Order2()


%% THIS HAS NEVER BEEN USED FOR THE PHD THESIS AND IS COPIED FROM "COMPUATATIONAL METHODS IN BIOMECHANICS" ASSINGMENT 3 %%% 
  clc;
  close all; 
  clear all;
  
  %------------------------------------------------------------------------
  % SETTINGS
  %------------------------------------------------------------------------
  % number of elements
  n_elem=2048;
  % number of grid points
  n=n_elem+1;
  % total number of grid points
  num_of_points = n;
  
  % Length
  L=1; %cm
  % grid spacing
  dx = L/n; % cm
  
  % end time
  t_end = 10.0; % ms
  % time step for the PDE
  time_step_pde = 0.01;
  % number of time steps for dynamic PDE solver
  num_of_steps_pde = t_end/time_step_pde;
  % number of ODE steps per one time step of PDE
  num_of_steps_ode=1;
  
  method='CN';
 
  %frequency of stimulation
  %f=60;
  
  % stimulation points
  %stim_pnts = zeros(n);
  % here, only one node in the middle of the domain is stimulated!!!
  stim_pnts(1) = ceil(num_of_points/2);
  %arbitrary point of stimulation
  %x_stim=0.1;
  %stim_pnts(1) = ceil(x_stim/dx);
  
  % start time of stimulation
  t_start_stimulation = 0.0;
  % stop time of stimulation
  t_end_stimulation = 0.5;
  
  %fast-twitch
  if(n_elem>16)
      I_stim=75*(n_elem/16); %the reference case must be adjusted by experimenting
  else
      I_stim=75;
  end
  
  %output time
  t_out=3;%ms  
  
  %------------------------------------------------------------------------
  %MATERIAL PARAMETERS
  %------------------------------------------------------------------------
   % surface area to volume ratio
  A_m = 500.0; % /cm
  % membrane capacitance, fast-twitch
  C_m = 1.0; % microF/cm^2
  % effective conductivity
  sigma_eff = 3.828; % mS/cm 
  % Diffusion coeff.
  D=sigma_eff/A_m/C_m;
  
  %------------------------------------------------------------------------
  % INITIALISE VARIABLES
  %------------------------------------------------------------------------
  % the time step number at which the stimulation starts
  if(t_start_stimulation ~= 0.0) 
    start_stim = t_start_stimulation/time_step_pde
  else
    start_stim=1;
  end
  % the time step number at which the stimulation stops
  stop_stim = t_end_stimulation/time_step_pde;
  
  % stimulation curtent
  i_Stim = zeros(num_of_points,num_of_steps_pde); 
  
  % membrane voltage as derived from the Hodgkin-Huxley model -- V*
  vS_hh = zeros(num_of_points,1);
  % transmembrane voltage after parabolic PDE evaluation -- V^{k+1}
  V_m = -75*ones(num_of_points,1);
  % transmembrane voltage at the output node for each time step
  V_m_time = zeros(num_of_steps_pde,num_of_points);
  
  %------------------------------------------------------------------------
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
  
  % Initialise constants and state variables for cellular model
  [CONSTANTS] = initConsts('HODGKIN_HUXLEY');
  [INIT_STATES] = initStates('HODGKIN_HUXLEY'); 
  
  % store all state variables for all points
  ALL_STATES = zeros(num_of_points,4); % 4 is the number of STATE variables
  % initialise ALL_STATES
  for i = 1:num_of_points
    ALL_STATES(i,:) = INIT_STATES;
  end
  
  %------------------------------------------------------------------------
  % STIMULATION
  %------------------------------------------------------------------------
  % set the stimulation at the specified nodes and times
  for i=start_stim:stop_stim
    for j=1:length(stim_pnts)
      node = stim_pnts(j);
      length(stim_pnts);
      i_Stim(node,i) = I_stim;
    end
  end 
  
  %------------------------------------------------------------------------
  %DISCRETIZATION OF PDE
  %------------------------------------------------------------------------
  % STIFFNESS MATRIX
  % first order dynamic problem -->  KK * V_m = bb
  KK= StiffnessMatrix(method,num_of_points,time_step_pde,dx,D); 
  
  %------------------------------------------------------------------------
  % SOLUTION PROCESS
  %------------------------------------------------------------------------ 
  % loop over the PDE time steps
  for time = 1:num_of_steps_pde   
    
    fprintf('step #   : %d \t', time);
    fprintf('time [ms]: %f \n', (time-1)*time_step_pde);
    
    % Set half of the timespan to integrate over the ODE
    tspan = (time-1):0.5/ (num_of_steps_ode):(time-0.5);
    tspan=tspan*time_step_pde;
    
    % Integrate the cellular model at each discretisation point
    [V_m,ALL_STATES] =SolveCellular(2,time,num_of_points,ALL_STATES,CONSTANTS,i_Stim,tspan); 
    
    %----------------------------------------------------------------------
    % PARABOLIC EQUATION
    bb=setRHS(method,V_m,time_step_pde,dx,D);
    % SOLVE THE PARABOLIC PDE
    V_m = KK\bb;
    
    %----------------------------------------------------------------------
    % update the cell models transmembrane voltage
    ALL_STATES(:,1)  = V_m; 
    
    % Set halsf of the timespan to integrate over the ODE
    tspan = (time-0.5):0.5/ num_of_steps_ode:time;
    tspan=tspan*time_step_pde;
    % Integrate the cellular model at each discretisation point
    [V_m,ALL_STATES] =SolveCellular(2,time,num_of_points,ALL_STATES,CONSTANTS,i_Stim,tspan); 
    
    %----------------------------------------------------------------------
    % update the cell models transmembrane voltage
    ALL_STATES(:,1)  = V_m;
    % store the transmembrane voltage at the output node
    V_m_time(time,:) = V_m;

  end %time
  
  OutToFile(V_m_time,t_out,time_step_pde,method); 
  
  tt = linspace(0,t_end,num_of_steps_pde);
  x_a = linspace(0,L,num_of_points);
  surf(x_a,tt,V_m_time);
  
end

function []=OutToFile(V_m_time,t_out,time_step,method)
time=t_out/time_step;

outfile=fopen(strcat('out','Vm_',num2str(time_step),'_',method,'.txt'),'w');
fprintf(outfile,'t_end=%f\ntime_step=%f\n',t_out,time_step);
fprintf(outfile,'%6.6f\n',V_m_time(time,:));
end
