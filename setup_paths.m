% setup_paths.m - Add all project directories to MATLAB/Octave path
%
% Run this script once before running any example:
%   >> setup_paths
%   >> test_multirate_simple

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'utils'));
addpath(fullfile(root, 'examples'));
fprintf('Paths added.\n');
