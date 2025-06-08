clear all
close all
clc

Y = struct( ...
    'N2', 0.767082, ...
    'O2', 0.232918, ...
    'NO', 2.69731e-16, ...
    'O',  1.25697e-41, ...
    'N',  2.39879e-80);

H01 = 18.59e6; % J/kg
 M1 = 6;

T1 = StaticTemperatureUpstream(H01, M1, Y);
disp('Static Temperature Upstream:');
disp(T1)
