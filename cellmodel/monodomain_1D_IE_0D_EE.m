%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This programme is intended to be a test case (convergence survey) for
% the time integrators (such as explicit Euler, Heun, implicit Euler etc.)
% in the framework of the 0D CellML Model, only. 
%
% The considered Problem is a system of coupled ODEs and given initial
% conditions, which needs to be integrated from t_start to t_end, both fix.
% The numerical integration can be done by different integration schemes
% (the 'time integrators' from above) as well as with variable time step
% sizes. We will then consider the numeric solution of the system at t_end
% for both, different schemes and time step sizes. 
%
% Within the survey, we will measure the effort of the schemes (CPU time)
% that they need to get to the numeric solution. Goal of the survey is to
% find a touple (scheme, time step size) which yields a numerical solution
% of the same quality (same order of the error) with less computational
% effort (shorter CPU time).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. set up the computation of 1 cell.
% 2. set up time measurement.
% 3. set up FOR loops and respective information vectors.
% 4. set up figures OR tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = monodomain_1D_IE_0D_EE()


%% THIS HAS NEVER BEEN USED FOR THE PHD THESIS AND IS COPIED FROM "COMPUATATIONAL METHODS IN BIOMECHANICS" ASSINGMENT 3 %%% 
  clc;
  close all; 
  clear all;
  
  format long g
  
  % Time step for PDE
  PDEtime_step = 0.025; % default: 0.1 ms change this?
  %------------------------------------------------------------------

  %------------------------------------------------------------------
  % STIMULATION
  Stim=1; % 0 = no, 1 = yes.
  if(Stim)
    stim_current=2000;
  else
    stim_current=0;
  end
  % initialize irregular vector with number of steps [1,2,4,8, ... ,32,50,(64,) ... 2^(max-2), 2^(max+3)]
  normal=0; % whether to use euler as best solution or heun
  Wied=1000;
  step_iter_max = 7; % mach 7
  SAVETHIS=step_iter_max;
  N_steps=zeros(1,step_iter_max);
  for step_iter=1:step_iter_max-1
    N_steps(step_iter)=2^(step_iter-1);
  end
  N_steps(end)=2^(step_iter_max+3); % was 6 for first (normal) execution
  if(step_iter_max>6)
      N_steps=[N_steps(1:6),50,N_steps(7:end)];
      step_iter_max=step_iter_max+1;
  end
  %------------------------------------------------------------------------
  
  % Initialise constants and state variables for cellular model
  [CONSTANTS0]   = initConsts('HODGKIN_HUXLEY');
  [INIT_STATES0] = initStates('HODGKIN_HUXLEY');
  %------------------------------------------------------------
  
  % amount of schemes that will be testesd:
  schemeCount = 3;
  %----------------------------------------
  
  % set up matrices for measurements..
  % A: state variable 1, B: state variable 2, C: state cariable 3, D: state
  % variable 4,  all of them at time PDEtime_step.
  % E: elapsed CPU time. Do post measurement for error stuff, so later 
  % there will be another matrix
  A=zeros(schemeCount,step_iter_max); B=A; C=A; D=A; E=A;
  
  for scheme_iter=1:schemeCount
      if(normal)
        % CARE: N_steps(end) will only be used for euler explicit.:-------
        if (scheme_iter==2) step_iter_max=step_iter_max-1; end
      else
        if (scheme_iter==1)
          step_iter_max=step_iter_max-1;
        elseif(scheme_iter==2)
          step_iter_max=step_iter_max+1;
        elseif(scheme_iter==3)
          step_iter_max=step_iter_max-1;
        end
      end
      %----------------------------------------------------------------- 
      
      for step_iter=1:step_iter_max
          N_GP = N_steps(step_iter)+1; % the number of grid points
          
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
  % algebraicVariableCount=10;

  % Set numerical accuracy options for ODE solver
  % options = odeset('RelTol', 1e-07, 'AbsTol', 1e-07, 'MaxStep', 0.1);  %           ??
 
  %---------------------------------------------------------------------
  % INTEGRATE CELLULAR MODEL
  
  % Set timespan to integrate over (equidistant time discretization)
          tspan = linspace(0,PDEtime_step,N_GP);
  % [VOI2, STATES2] = ode15s(@(VOI, STATES)computeRates_HODGKIN_HUXLEY(VOI, STATES, CONSTANTS, stim_current), tspan2, INIT_STATES, options);
          if(scheme_iter==1)
              tic;
              for high=1:Wied
                  STATES = ode1(@(VOI, STATES)computeRates_HODGKIN_HUXLEY(VOI, STATES, CONSTANTS0, stim_current), tspan, INIT_STATES0);
              end
              el = toc;
          elseif(scheme_iter==2)
              tic;
              for high=1:Wied
                  STATES = ode2(@(VOI, STATES)computeRates_HODGKIN_HUXLEY(VOI, STATES, CONSTANTS0, stim_current), tspan, INIT_STATES0);
              end
              el = toc;
          elseif(scheme_iter==3)
              tic;
              for high=1:Wied
                  STATES = ode3(@(VOI, STATES)computeRates_HODGKIN_HUXLEY(VOI, STATES, CONSTANTS0, stim_current), tspan, INIT_STATES0);
              end
              el = toc;
          else
              tic;
              for high=1:Wied
                  STATES = ode4(@(VOI, STATES)computeRates_HODGKIN_HUXLEY(VOI, STATES, CONSTANTS0, stim_current), tspan, INIT_STATES0);
              end
              el = toc;
          end
          
          A(scheme_iter,step_iter) = STATES(end,1);
          B(scheme_iter,step_iter) = STATES(end,2);
          C(scheme_iter,step_iter) = STATES(end,3);
          D(scheme_iter,step_iter) = STATES(end,4);
          E(scheme_iter,step_iter) = el;
%? Compute algebraic variables
%?      [RATES, ALGEBRAIC] = computeRates_HODGKIN_HUXLEY(VOI, STATES, CONSTANTS, i_stim(point,time));
%?      ALGEBRAIC = computeAlgebraic_HODGKIN_HUXLEY(ALGEBRAIC, CONSTANTS, STATES, VOI, I_HH(point,time));


      end
  end
  %tt = linspace(0,t_end,num_of_dt);
  %x_a = linspace(0,1,num_of_points);
  %surf(x_a,tt,V_m_time);
  
  % post processing: ------------------------------------------------------
  % compute relative errors (assume the false assumption that M(1,end)
  % holds the most accurate solution (thats the very fine discretization
  % with euler explicit, ode1))
  
  ErrA=zeros(schemeCount,SAVETHIS); ErrB=ErrA; ErrC=ErrA; ErrD=ErrA;
  if(normal)
    ErrA(:,:)=abs((A(:,1:end-1)-A(1,end))./A(1,end));
    ErrB(:,:)=abs((B(:,1:end-1)-B(1,end))./B(1,end));
    ErrC(:,:)=abs((C(:,1:end-1)-C(1,end))./C(1,end));
    ErrD(:,:)=abs((D(:,1:end-1)-D(1,end))./D(1,end));
  else
    ErrA(:,:)=abs((A(:,1:end-1)-A(2,end))./A(2,end));
    ErrB(:,:)=abs((B(:,1:end-1)-B(2,end))./B(2,end));
    ErrC(:,:)=abs((C(:,1:end-1)-C(2,end))./C(2,end));
    ErrD(:,:)=abs((D(:,1:end-1)-D(2,end))./D(2,end));
  end
  %------------------------------------------------------------------------
  
  % short output
  E, ErrA, ErrB, ErrC, ErrD
  save('StimEl_relERR_A_B_C_Dunnormal.mat','E','ErrA','ErrB','ErrC','ErrD');
  save('StimEl_relERR_A_B_C_Dunnormal.txt','E','ErrA','ErrB','ErrC','ErrD','-ascii');
  
end

function []=OutToFile(V_m_time,t_out,time_step,method)
time=t_out/time_step;

outfile=fopen(strcat('out','Vm_',num2str(time_step),'_',method,'.txt'),'w');
fprintf(outfile,'t_end=%f\ntime_step=%f\n',t_out,time_step);
fprintf(outfile,'%6.6f\n',V_m_time(time,:));
end
