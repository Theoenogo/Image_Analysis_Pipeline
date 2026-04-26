% masterScript.m

% Define mainFolderPath as a global variable and set its value
global mainFolderPath;
msgbox('Please select the folder containing the main path.', 'Prompt', 'modal');
mainFolderPath = uigetdir('Select the main folder path');
if isequal(mainFolderPath, 0) % User pressed cancel
    disp('Folder selection was canceled. Exiting script.');
    return;
end

% Run the first script
A_renameFolders_v2();

% Run the second script
B_deleteScanProtocolFiles();

% Run the third script
C2_green_stack_offset_Cy5_recursive();

% Run the fourth script
D2_cy5_stack_no_offset_recursive();

% Run the fifth script
E_move_decon_folders();

% Gather user input for runDeconvolutionJobGFP
msgbox('Please select the PSF image for GFP.', 'Prompt', 'modal');
[psfFileGFP, psfPathGFP] = uigetfile('*.tif', 'Select the PSF image for GFP');
psfImagePathGFP = fullfile(psfPathGFP, psfFileGFP);

% Gather user input for runDeconvolutionJob
msgbox('Please select the PSF image for Cy.', 'Prompt', 'modal');
[psfFileCy, psfPathCy] = uigetfile('*.tif', 'Select the PSF image for Cy');
psfImagePathCy = fullfile(psfPathCy, psfFileCy);
msgbox('Please select the Decon folder.', 'Prompt', 'modal');
baseFolderPath = uigetdir('', 'Select the base folder for both GFP and Cy (e.g., Decon)');

% Detect the operating system
if ismac
    osType = 'MACI64';
elseif isunix
    osType = 'UNIX'; % Add more specific checks if needed
elseif ispc
    osType = 'PCWIN'; % Add more specific checks if needed
else
    error('Unsupported OS');
end

disp(['Detected OS: ', osType]); % Add this line for debugging

% Update paths for deconvolution scripts
deconvolutionScriptFolder = '/Users/theodorecarter/Documents/Data_Processing/Key_Scripts/MatLab/Recursive_Streamlined';

% Construct and execute the command based on the OS
if strcmpi(osType, 'PCWIN')  % Windows
    command1 = sprintf('start matlab -r "addpath(''%s''); F_GFP_runDeconvolutionJob_parallel(''%s'', ''%s'')"', deconvolutionScriptFolder, psfImagePathGFP, baseFolderPath);
    command2 = sprintf('start matlab -r "addpath(''%s''); F_Cy5_runDeconvolutionLab2_parallel(''%s'', ''%s'')"', deconvolutionScriptFolder, psfImagePathCy, baseFolderPath);
elseif strcmpi(osType, 'MACI64')  % macOS
    command1 = sprintf('open -n /Applications/MATLAB_R2024a.app --args -r "addpath(''%s''); F_GFP_runDeconvolutionJob_parallel(''%s'', ''%s'')"', deconvolutionScriptFolder, psfImagePathGFP, baseFolderPath);
    command2 = sprintf('open -n /Applications/MATLAB_R2024a.app --args -r "addpath(''%s''); F_Cy5_runDeconvolutionLab2_parallel(''%s'', ''%s'')"', deconvolutionScriptFolder, psfImagePathCy, baseFolderPath);
else
    error('Unsupported OS: %s', osType);
end

system(command1);
system(command2);

% Display message
disp('Scripts are being executed in separate MATLAB instances.');

% Pause to prevent system sleep during deconvolution (after deconvolution commands)
pause(60); % Pause for 60 seconds (adjust as needed)