%% Assignment 4 - Ferranti Effect
% 220 kV, 50 Hz, 200 km transmission line
% Analytical nominal-pi benchmark for Simulink verification

clear;
clc;
close all;

fprintf('============================================================\n');
fprintf(' Assignment 4: Ferranti Effect - 220 kV, 50 Hz, 200 km\n');
fprintf('============================================================\n\n');

%% INPUT DATA

VLL_kV = input('Line-line voltage in kV [220]: ');
if isempty(VLL_kV)
    VLL_kV = 220;
end

f = input('Frequency in Hz [50]: ');
if isempty(f)
    f = 50;
end

len_km = input('Line length in km [200]: ');
if isempty(len_km)
    len_km = 200;
end

R_km = input('Series resistance in ohm/km/phase [0.05]: ');
if isempty(R_km)
    R_km = 0.05;
end

L_km_mH = input('Series inductance in mH/km/phase [1]: ');
if isempty(L_km_mH)
    L_km_mH = 1;
end

C_km_uF = input('Shunt capacitance in uF/km/phase [0.01]: ');
if isempty(C_km_uF)
    C_km_uF = 0.01;
end

P_light_MW = input('Light-load active power in MW [50]: ');
if isempty(P_light_MW)
    P_light_MW = 50;
end

pf_light = input('Light-load lagging power factor [0.95]: ');
if isempty(pf_light)
    pf_light = 0.95;
end

%% CONVERT LINE PARAMETERS

VLL  = VLL_kV * 1e3;
w    = 2 * pi * f;

R = R_km * len_km;
L = (L_km_mH * 1e-3) * len_km;
C = (C_km_uF * 1e-6) * len_km;

Vphs = VLL / sqrt(3);
Vs   = Vphs;

%% NOMINAL-PI ABCD PARAMETERS

Z    = R + 1j * w * L;
Y    = 1j * w * C;

A    = 1 + Z * Y / 2;
B    = Z * (1 + Z * Y / 4);
Cabc = Y * (1 + Z * Y / 4);
D    = A;

%% NO-LOAD CONDITION

Vr0 = Vs / A;

Ir0 = Cabc * Vr0;
Is0 = A * Ir0 + Cabc * Vr0;
S0  = 3 * Vs * conj(Is0);

Vrecv0_kV = abs(Vr0) * sqrt(3) / 1e3;

Ferranti_percent = ...
    (Vrecv0_kV - VLL_kV) / VLL_kV * 100;

%% LIGHT-LOAD CONDITION

P = P_light_MW * 1e6;
Q = P * tan(acos(pf_light));

Vr = Vphs;

for k = 1:1000

    Ir = conj((P + 1j * Q) / (3 * Vr));
    Vr_new = (Vs - B * Ir) / A;

    if abs(Vr_new - Vr) < 1e-8
        Vr = Vr_new;
        break;
    end

    Vr = Vr_new;

end

Ir_light = conj((P + 1j * Q) / (3 * Vr));

Is_light = A * Ir_light + Cabc * Vr;
SL       = 3 * Vs * conj(Is_light);

VrecvL_kV = abs(Vr) * sqrt(3) / 1e3;

Light_deviation_percent = ...
    (VrecvL_kV - VLL_kV) / VLL_kV * 100;

%% SHUNT-REACTOR MITIGATION

% Tune the receiving-end reactor so that the no-load
% receiving voltage is restored to the nominal 220 kV.

target = VLL;

fun = @(Br) abs( ...
    Vs / (A + B * (-1j * Br)) ...
    ) * sqrt(3) - target;

Br = fzero(fun, [0 0.005]);

Yr  = -1j * Br;
VrR = Vs / (A + B * Yr);

IrR = Yr * VrR;
IsR = A * IrR + Cabc * VrR;
SR  = 3 * Vs * conj(IsR);

Qreactor_Mvar = ...
    3 * (VLL / sqrt(3))^2 * Br / 1e6;

Lreactor_H = 1 / (w * Br);

VrecvR_kV = abs(VrR) * sqrt(3) / 1e3;

Reactor_deviation_percent = ...
    (VrecvR_kV - VLL_kV) / VLL_kV * 100;

%% DISPLAY LINE PARAMETERS

fprintf('\n--- LINE PARAMETERS ---\n');

fprintf('R = %.4f ohm/phase\n', R);
fprintf('L = %.4f H/phase\n', L);
fprintf('C = %.6e F/phase\n', C);
fprintf('Z = %.4f + j%.4f ohm/phase\n', ...
    real(Z), imag(Z));
fprintf('Y = j%.8f S/phase\n', imag(Y));

%% DISPLAY NO-LOAD RESULTS

fprintf('\n--- NO LOAD ---\n');

fprintf('Receiving V = %.3f kV L-L\n', Vrecv0_kV);
fprintf('Ferranti voltage rise = %.3f %%\n', ...
    Ferranti_percent);
fprintf('Receiving current = %.3f A\n', abs(Ir0));

%% DISPLAY LIGHT-LOAD RESULTS

fprintf('\n--- LIGHT LOAD ---\n');

fprintf('P = %.2f MW, pf = %.3f lagging\n', ...
    P_light_MW, pf_light);

fprintf('Qload = %.3f MVAr\n', Q / 1e6);
fprintf('Receiving V = %.3f kV L-L\n', VrecvL_kV);
fprintf('Deviation = %.3f %%\n', Light_deviation_percent);

%% DISPLAY MITIGATION RESULTS

fprintf('\n--- MITIGATION ---\n');

fprintf('Reactor rating = %.3f MVAr at 220 kV\n', ...
    Qreactor_Mvar);

fprintf('Equivalent reactor L = %.4f H/phase\n', ...
    Lreactor_H);

fprintf('Compensated receiving V = %.3f kV L-L\n', ...
    VrecvR_kV);

%% RESULTS TABLE

Case = {
    'No load'
    'Light load'
    'No load + reactor'
    };

Vrecv_kV = [
    Vrecv0_kV
    VrecvL_kV
    VrecvR_kV
    ];

Deviation_pct = [
    Ferranti_percent
    Light_deviation_percent
    Reactor_deviation_percent
    ];

disp(table(Case, Vrecv_kV, Deviation_pct));

%% RECEIVING VOLTAGE VS LOAD

Pscan_MW = 0:5:100;
Vscan_kV = zeros(size(Pscan_MW));

for n = 1:numel(Pscan_MW)

    if Pscan_MW(n) == 0

        Vscan_kV(n) = Vrecv0_kV;

    else

        Pn  = Pscan_MW(n) * 1e6;
        Qn  = Pn * tan(acos(pf_light));
        Vrn = Vphs;

        for k = 1:1000

            In = conj((Pn + 1j * Qn) / (3 * Vrn));
            Vnew = (Vs - B * In) / A;

            if abs(Vnew - Vrn) < 1e-8
                Vrn = Vnew;
                break;
            end

            Vrn = Vnew;

        end

        Vscan_kV(n) = abs(Vrn) * sqrt(3) / 1e3;

    end

end

%% PLOT 1: VOLTAGE VS LOAD

figure('Color', 'w');

plot(Pscan_MW, Vscan_kV, ...
    'LineWidth', 1.8);

hold on;

yline(VLL_kV, '--', '220 kV');

grid on;

xlabel('Receiving-end active power (MW)');
ylabel('Receiving-end voltage (kV L-L)');

title('Ferranti Effect: Receiving-End Voltage vs Load');

legend( ...
    'Calculated voltage', ...
    'Nominal voltage', ...
    'Location', 'best');

%% PLOT 2: MITIGATION COMPARISON

figure('Color', 'w');

bar(categorical(Case), Vrecv_kV);

grid on;

ylabel('Receiving-end voltage (kV L-L)');

title('Ferranti Effect and Shunt-Reactor Mitigation');
