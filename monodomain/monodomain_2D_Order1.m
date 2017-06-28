%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = monodomain_1D_Order1()


%% THIS HAS NEVER BEEN USED FOR THE PHD THESIS AND IS COPIED FROM "COMPUATATIONAL METHODS IN BIOMECHANICS" ASSINGMENT 3 %%% 
  clc;
  close all; 
  clear all;
  
  %------------------------------------------------------------------------
  % SETTINGS
  %------------------------------------------------------------------------
  %number of elements in each direction
  n_elem=24;
  % number of grid points in each direction
  n=n_elem+1;
  % total number of grid points
  num_of_points = n*n;
  
  % Length
  L=1; %cm
  % grid spacing
  dx = L/n; % cm
  
   % stop time
  t_end = 10.0; % ms
  % time step for the PDE
  time_step_pde = 0.05;
  % number of time steps for dynamic PDE solver
  num_of_steps_pde = t_end/time_step_pde;
  % number of ODE runs per one time step of PDE
  num_of_steps_ode=50;
  
  method='ImplicitEuler';
 
  %frequency of stimulation
  %f=60;
  
  % stimulation points
  %stim_pnts = zeros(n);
  % here, only one node in the middle of the domain is stimulated!!!
  stim_pnts(1) = ceil(num_of_points/2)%check if the point s not on the boundary
  %stim_pnts(1) = 25;
  %x_stim=0.1;
  %stim_pnts(1) = ceil(x_stim/dx);
  
  % start time of stimulation
  t_start_stimulation = 0.0;
  % stop time of stimulation
  t_end_stimulation = 0.1;
  
  %fast twitch
  I_stim=2000*(n/24)^2;
  
  %output time
  t_out=3;%ms
  
  %output node
  out_node=13;
  
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
    start_stim = t_start_stimulation/time_step_pde;
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
  % V_m at the output node
  V_m_out= zeros(num_of_steps_pde,1);
  figure(1);
  writerObj=VideoWriter('Monodomain.avi');
  open(writerObj)
  
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
  KK= StiffnessMatrix_2D(method,num_of_points,time_step_pde,dx,D); 
  
  %---------------------------------------------------------------------
  % SOLUTION PROCESS
  %--------------------------------------------------------------------- 
  % loop over the PDE time steps
  for time = 1:num_of_steps_pde
    
    %fprintf('step #   : %d \t', time);
    %fprintf('time [ms]: %f \n', (time-1)*time_step);
    
    % Set timespan to integrate over the ODE
    tspan = (time-1):1/ num_of_steps_ode:time;
    tspan=tspan*time_step_pde;
    
    % Integrate the cellular model at each discretisation point
    for point = 1:num_of_points
      
      % load the last state as initial state for the next integration
      LAST_STATES = ALL_STATES(point,:);
      
      %---------------------------------------------------------------------
      % INTEGRATE CELLULAR MODEL WITH ODE SOLVER     
      [STATES] = ode1(@(VOI, STATES)computeRates_HODGKIN_HUXLEY(VOI, STATES, CONSTANTS, i_Stim(point,time)), tspan, LAST_STATES);

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
    
    % PARABOLIC EQUATION
    bb=setRHS_2D(method,V_m,time_step_pde,dx,D); 
    % SOLVE THE PARABOLIC PDE
    V_m = KK\bb;
    
    %---------------------------------------------------------------------
    % update the cell models transmembrane voltage
    ALL_STATES(:,1)  = V_m;
    %store the transmembrane voltage for all time steps
    V_m_time(time,:) = V_m;
    % store the transmembrane voltage at the output node
    V_m_out(time)=V_m(out_node);
    
    y_a = (0:1/(n-1):1);
    x_a = y_a;
    u = zeros(n,n);
    for i = 1:n
      u(:,i) = V_m((i-1)*n+1:i*n);
    end
    %surf(meshgrid(x_a)',meshgrid(y_a),u);
    %figure(1);
    clf;
    surf(meshgrid(x_a)',meshgrid(y_a),u);
    axis([0 1 0 1 -100 40 0 1]);
    %refresh(1);
    drawnow;
    %MV(time)=getframe;
    thisimage=getframe;
    writeVideo(writerObj,thisimage);

   end %time
  
  OutToFile(V_m_time,t_out,time_step_pde,method); 
  %figure
  %movie(MV,num_of_steps_pde);
  
  tt = linspace(0,t_end,num_of_steps_pde);
%   x_a = linspace(0,L,n);
%   y_a=linspace(0,L,n);
%   u = zeros(n,n);
%   for i = 1:n
%       u(:,i) = V_m_time(t_end,(i-1)*n+1:i*n);
%   end
% figure(time);
% surf(meshgrid(x_a)',meshgrid(y_a),u);
close(writerObj);
V_m_out  
figure(999);
plot(tt, V_m_out);
  
end

function []=OutToFile(V_m_time,t_out,time_step,method)
time=t_out/time_step;

outfile=fopen(strcat('out','Vm_',num2str(time_step),'_',method,'.txt'),'w');
fprintf(outfile,'t_end=%f\ntime_step=%f\n',t_out,time_step);
fprintf(outfile,'%6.6f\n',V_m_time(time,:));
end

%function []=Movie
