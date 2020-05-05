%% Compute ALE for Chronic Pain Meta-analysis Project
% by AX 11/15/2019
% -------------------------------------
% This script is a higher-level script that runs all operations pertaining to
% computing ALE (using Eickhoff's scripts)
%   a) thresholds the main effects clusters at z>3.09
%   b) generates the negative & positive contrast images for contrasts
%   c) thresholds the contrast clusters at a liberal threshold of z>1.6 and k>50
%       % N.B. : this will eventually get re-thresholded with FSL at z>3.09
%   d) thresholds the conjunction clusters at z>3.09 and k>15
% -------------------------------------

%% Setup

curPath=genpath(pwd); addpath(curPath); % add the current path
cd ale % go to ALE folder

%% Main and Sub-analyses 01/13/2020

% Input
ale_inputCoords('aberrant_20200113.xls'); % Aberrant
ale_inputCoords('painCoords_20200113.xls'); % Patients (signed)

% Compute
ale_estimateALE('pain_20200113.xlsx');

%% Subset Perceptual Pain 01/21/2020

% Input
ale_inputCoords('perceptualData_2020-01-21.xls');

% Compute (not Control > Patient because n=8)
ale_estimateALE('perceptual_20200121.xlsx');
