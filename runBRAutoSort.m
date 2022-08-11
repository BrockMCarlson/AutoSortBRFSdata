% runBRAutoSort.m
%
%


clear
close all
PostSetup('BrockHome')

DATADIR = 'D:\all BRFS\151222_E';
cd(DATADIR);

filename = 'D:\all BRFS\151222_E\151222_E_brfs001';
multiplier = 2.2;
[ppNEV, WAVES] = offlineBRAutoSort(filename,multiplier);
