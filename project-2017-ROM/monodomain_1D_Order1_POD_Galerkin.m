%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = monodomain_1D_Order1_POD_Galerkin()


%% THIS HAS NEVER BEEN USED FOR THE PHD THESIS AND IS COPIED FROM "COMPUATATIONAL METHODS IN BIOMECHANICS" ASSINGMENT 3 %%% 
  clc;
  close all; 
  clear all;
  
  %------------------------------------------------------------------------
  % SETTINGS
  %------------------------------------------------------------------------
  %number of elements
  n_elem=16;
  % number of grid points
  n=n_elem+1;
  % total number of grid points
  num_of_points = n;
  
  % Length
  L=1; %cm
  % grid spacing
  dx = L/n; % cm
  
   % end time
  t_end = 5.0; % ms
  % time step for the PDE
  time_step_pde = 0.01;
  % number of time steps for dynamic PDE solver
  num_of_steps_pde = t_end/time_step_pde;
  % number of ODE runs per one time step of PDE
  num_of_steps_ode=1;
  
  method='ImplicitEuler';
 
  %frequency of stimulation
  %f=60;
  
  % stimulation points
  %stim_pnts = zeros(n);
  % here, only one node in the middle of the domain is stimulated!!!
  stim_pnts(1) = ceil(num_of_points/2);
  % arbitrary stimulation point
  %x_stim=0.1;
  %stim_pnts(1) = ceil(x_stim/dx);
  
  % start time of stimulation
  t_start_stimulation = 0.0;
  % stop time of stimulation
  t_end_stimulation = 0.5;
  
  %fast twitch
  if(n_elem>16)
      I_stim=75*(n_elem/16); %the reference case must be adjusted by experimenting
  else
      I_stim=75;
  end
  %output time
  t_out=3;%ms 
  
  % arbitrary output nodes
  out_node=ceil(num_of_points/2);
  
  % arbitrary output nodes to measure the propagation velocity
  L1=L/8;
  dL=L/4;
  L2=L1+dL;
  out_node_1=round(L1/L*n_elem)
  out_node_2=round(L2/L*n_elem)+1
  
  %------------------------------------------------------------------------
  %MATERIAL PARAMETERS
  %------------------------------------------------------------------------
   % surface area to volume ratio
  A_m = 500.0; % /cm
  % membrane capacitance, fast-twitch
  C_m = 1.0; % microF/cm^2
  % effective conductivity
  sigma_eff = 3.828; %3.828; % mS/cm 
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
  % R elem n * t
  i_Stim = zeros(num_of_points,num_of_steps_pde); 
  
  % membrane voltage as derived from the Hodgkin-Huxley model -- V*
  vS_hh = zeros(num_of_points,1);
  % transmembrane voltage after parabolic PDE evaluation -- V^{k+1}
  % y0
  % R elem n
  V_m = -75*ones(num_of_points,1);
  % Uk laden
  % Basisvektoren POD
  % R elem n * k
  Uk = csvread('Ukt.csv');
  [m,n] = size(Uk);
  % R elem k
  %Anzahl Basisvektoren k
  num_of_points_tilde = n;
  % y~0 = Ukt * y0
  % R elem k
  %V_m_tilde = Uk' * V_m
  % transmembrane voltage at the output node for each time step
  % R elem t * n
  V_m_time = zeros(num_of_steps_pde,num_of_points);
  % R elem t * k
  V_m_time_tilde = zeros(num_of_steps_pde,num_of_points_tilde);
  
  
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
  %------------------------------------------------------------------
  % STIFFNESS MATRIX
  % first order dynamic problem -->  KK * V_m = bb
  % R elem n * n
  % TODO: Boundary Conditions?
  %KK= StiffnessMatrix_POD_Galerkin(method,num_of_points,time_step_pde,dx,D);
  KK= StiffnessMatrix(method,num_of_points,time_step_pde,dx,D);
  % R elem k * k
  KK_tilde = Uk' * KK * Uk;
  
  % Von Neumann hier verwendet. Andere M�glichkeit: in KK (auskommentiert in Methode Stiffness Matrix)
  % homogeneous Neumann B.C. - no flow
  % right side
    %KK_tilde(num_of_points_tilde,num_of_points_tilde) = 1;
    %KK_tilde(num_of_points_tilde,num_of_points_tilde-1) = -1;
  % left side
    %KK_tilde(1,1) = 1;
    %KK_tilde(1,2) = -1; 
  
  %---------------------------------------------------------------------
  % SOLUTION PROCESS
  %--------------------------------------------------------------------- 
  % loop over the PDE time steps
  for time = 1:num_of_steps_pde
    
    fprintf('step #   : %d \t', time);
    fprintf('time [ms]: %f \n', (time-1)*time_step_pde);
    
    % Set timespan to integrate over the ODE
    tspan = (time-1):1/ num_of_steps_ode:time;
    tspan=tspan*time_step_pde;
    
    % Integrate the cellular model at each discretisation point
    % 
    [V_m,ALL_STATES] =SolveCellular(1,time,num_of_points,ALL_STATES,CONSTANTS,i_Stim,tspan);
    
    % aktualisiere V_m_tilde
    % Ukt * F(Uk*y(t)_tilde)
    %V_m_tilde = Uk' * V_m;
    
    % PARABOLIC EQUATION
    % Achtung: funktioniert bisher nur mit 1. Ordnung
    % R elem 
    % mit Grenzbedingungen
    %bb=setRHS(method,V_m_tilde,time_step_pde,dx,D);
    bb=setRHS(method,V_m,time_step_pde,dx,D);
    bb_tilde =Uk' * bb;
    % SOLVE THE PARABOLIC PDE
    % kk_tilde * Vm^k+1_tilde = Vm*_tilde nach Vm^k+1_tilde aufl�sen
    %  (R elem k * k) * ( R elem k) = ( R elem k) 
    V_m_tilde = KK_tilde\bb_tilde;
    
    %---------------------------------------------------------------------
    % update the cell models transmembrane voltage
    % F(Uk* y(t)_tilde)
    %R elem n = (R elem n * k) * (R elem k)
    ALL_STATES(:,1)  = Uk * V_m_tilde;
    % store the transmembrane voltage at the output node
    V_m_time(time,:) = Uk * V_m_tilde;
    % Rohdaten
    V_m_time_tilde(time, :) = V_m_tilde;

  end %time
  
  %OutToFile(V_m_time,  V_m_time_tilde, t_out,time_step_pde,method);
  OutToFile_time_space(V_m_time,time_step_pde,n_elem);
  OutNodeinTime(V_m_time,out_node_1,time_step_pde,n_elem);
  OutNodeinTime(V_m_time,out_node_2,time_step_pde,n_elem);
  
  tt = linspace(0,t_end,num_of_steps_pde);
  x_a = linspace(0,L,num_of_points);
  surf(x_a,tt,V_m_time);
  
end

function []=OutToFile(V_m_time, V_m_time_tilde, t_out,time_step,method)
time=t_out/time_step;
csvwrite('V_m_time_Galerkin.csv',V_m_time);
csvwrite('V_m_time_tilde_Galerkin.csv', V_m_time_tilde);
outfile=fopen(strcat('out','Vm_POD',num2str(time_step),'_',method,'.txt'),'w');
fprintf(outfile,'%6.6f\n',V_m_time(time,:));
end

function []=OutToFile_time_space(V_m_time,timestep,n_elem)
outfile=strcat('out','Vm_POD_time_space','_dt_',num2str(timestep),'_n_',num2str(n_elem),'_O1.csv');
csvwrite(outfile,V_m_time); 
end


function[]=OutNodeinTime(Vm_time,indx,timestep,n_elem)
outfile=fopen(strcat('out_Vm_POD_node_',num2str(indx),'_dt_',num2str(timestep),'_n_',num2str(n_elem),'_O1.txt'),'w');
fprintf(outfile,'%6.6f\n',Vm_time(:,indx));
fclose(outfile);
end

function []=PrintOutNodes(L,n_elem)
    outfile=fopen(strcat('out_L_',num2str(L),'_n_',num2str(n_elem),'.txt'),'w');
    dx=L/n_elem;
    x=0.0:dx:L;
    fprintf(outfile,'%6.6f\n',x);
    fclose(outfile);
end
