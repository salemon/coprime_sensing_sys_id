% setup_paths.m - Add all project directories to MATLAB/Octave path
%
% Run this script once before running any example:
%   >> setup_paths
%   >> validate_camera_ready_testn2

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'utils'));
addpath(fullfile(root, 'examples'));
fprintf('Paths added.\n');
