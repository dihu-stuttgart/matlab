
function [VOI, STATES, ALGEBRAIC, CONSTANTS] = mainFunction()
    % This is the "main function".  In Matlab, things work best if you rename this function to match the filename.
   [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel();
end

function [algebraicVariableCount] = getAlgebraicVariableCount()
    % Used later when setting a global variable with the number of algebraic variables.
    % Note: This is not the "main method".
    algebraicVariableCount =1;
end
% There are a total of 2 entries in each of the rate and state variable arrays.
% There are a total of 6 entries in the constant variable array.
%

function [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel()
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  algebraicVariableCount = getAlgebraicVariableCount();
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;

    % Set timespan to solve over
    tspan = [0, 10];

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 1);

    % Solve model with ODE solver
    [VOI, STATES] = ode15s(@(VOI, STATES)computeRates(VOI, STATES, CONSTANTS), tspan, INIT_STATES, options);

    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI);

    % Plot state variables against variable of integration
    [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends();
    figure();
    plot(VOI, STATES);
    xlabel(LEGEND_VOI);
    l = legend(LEGEND_STATES);
    set(l,'Interpreter','none');
end

function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_VOI = ''; LEGEND_CONSTANTS = '';
    LEGEND_VOI = strpad('t in component Aliev_Panfilov (dimensionless)');
    LEGEND_STATES(:,1) = strpad('V_m in component Aliev_Panfilov (dimensionless)');
    LEGEND_STATES(:,2) = strpad('r in component Aliev_Panfilov (dimensionless)');
    LEGEND_ALGEBRAIC(:,1) = strpad('I_HH in component Aliev_Panfilov (dimensionless)');
    LEGEND_CONSTANTS(:,1) = strpad('T in component Aliev_Panfilov (dimensionless)');
    LEGEND_CONSTANTS(:,2) = strpad('k in component Aliev_Panfilov (dimensionless)');
    LEGEND_CONSTANTS(:,3) = strpad('a in component Aliev_Panfilov (dimensionless)');
    LEGEND_CONSTANTS(:,4) = strpad('e_0 in component Aliev_Panfilov (dimensionless)');
    LEGEND_CONSTANTS(:,5) = strpad('m_1 in component Aliev_Panfilov (dimensionless)');
    LEGEND_CONSTANTS(:,6) = strpad('m_2 in component Aliev_Panfilov (dimensionless)');
    LEGEND_RATES(:,1) = strpad('d/dt V_m in component Aliev_Panfilov (dimensionless)');
    LEGEND_RATES(:,2) = strpad('d/dt r in component Aliev_Panfilov (dimensionless)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end

function [STATES, CONSTANTS] = initConsts()
    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];
    STATES(:,1) = 0;
    STATES(:,2) = 0;
    CONSTANTS(:,1) = 0.2;
    CONSTANTS(:,2) = 128;
    CONSTANTS(:,3) = 0.15;
    CONSTANTS(:,4) = 0.002;
    CONSTANTS(:,5) = 0.2;
    CONSTANTS(:,6) = 0.3;
    if (isempty(STATES)), warning('Initial values for states not set');, end
end

function [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS)
    global algebraicVariableCount;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        ALGEBRAIC = zeros(1, algebraicVariableCount);
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        ALGEBRAIC = zeros(statesRowCount, algebraicVariableCount);
        RATES = zeros(statesRowCount, statesColumnCount);
        utilOnes = ones(statesRowCount, 1);
    end
    RATES(:,2) =  (CONSTANTS(:,4)+( CONSTANTS(:,5).*STATES(:,2))./(CONSTANTS(:,6)+STATES(:,1))).*( - STATES(:,2) -  CONSTANTS(:,2).*STATES(:,1).*((STATES(:,1) - CONSTANTS(:,3)) - 1.00000));
    ALGEBRAIC(:,1) = piecewise({VOI>=10.0000&VOI<=10.5000, 10.0000 }, 0.00000);
    RATES(:,1) = (  - CONSTANTS(:,2).*STATES(:,1).*(STATES(:,1) - CONSTANTS(:,3)).*(STATES(:,1) - 1.00000) -  STATES(:,1).*STATES(:,2))+ALGEBRAIC(:,1);
   RATES = RATES';
end

% Calculate algebraic variables
function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI)
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        utilOnes = ones(statesRowCount, 1);
    end
    ALGEBRAIC(:,1) = piecewise({VOI>=10.0000&VOI<=10.5000, 10.0000 }, 0.00000);
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
