% There are a total of 4 entries in each of the rate and state variable arrays.
% There are a total of 8 entries in the constant variable array.

function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends_HODGKIN_HUXLEY()
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