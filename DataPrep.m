% The last dataset from Timon on 18 Nov 2020
clear; clear all; clc; close all;

% PG, VG points of data set
PGVG_table = readtable('mvnd_2e5_samples.csv');
PG_array = table2array(PGVG_table(:,[1:4]));
max_PG = max(PG_array) % maximum element of each column
min_PG = min(PG_array) 
VG_array = table2array(PGVG_table(:,[5:9]));
max_VG = max(VG_array) 
min_VG = min(VG_array) 

% stability (1) or not (0) of points
output_struct = load('reclassification.mat');
output_cell = struct2cell(output_struct);
output_array = cell2mat(output_cell);

% check PG points
PG_output = horzcat(PG_array,output_array);
output_mod = output_array;

% Print lines where Pg<0 (not allowed) and still classified as stable
count  = 0;
for i = 1:size(PG_array,1)
    for j = 1:size(PG_array,2)
        if PG_array(i,j)<0 && output_array(i) > 0
            PG_output(i,:);
            output_mod(i) = 0; % if any PG negative, then rewrite that point is not stable
            count = count + 1;
        end 
    end 
end
% The smallest PG line: [0.7311    0.4263    0.4522   -0.0037    1.0000] ->
% -3.7kW!

% Writing csv file
writematrix(table2array(PGVG_table),'NN_input.csv');
writematrix(output_mod,'NN_output_mod.csv');

