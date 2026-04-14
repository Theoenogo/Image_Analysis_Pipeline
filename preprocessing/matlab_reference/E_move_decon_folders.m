function E_move_decon_folders()

    % Access the global variable mainFolderPath
    global mainFolderPath;

    % Create a new folder named "Decon" in the main folder
    deconFolder = fullfile(mainFolderPath, 'Decon');
    mkdir(deconFolder);

    % Get a list of all subdirectories in the main folder
    subdirs = dir(mainFolderPath);
    subdirs = subdirs([subdirs.isdir]);  % Filter out non-directory entries
    subdirs = subdirs(~ismember({subdirs.name}, {'.', '..'}));  % Exclude "." and ".." entries

    % Loop through each subdirectory
    for i = 1:numel(subdirs)
        subdirName = subdirs(i).name;
    
    % Check if the subdirectory is named "Decon"
    if strcmp(subdirName, 'Decon')
        continue;  % Skip the "Decon" folder itself
    end
    
    subdirPath = fullfile(mainFolderPath, subdirName);
    
    % Check if the subdirectory contains a folder named "Decon"
    if exist(fullfile(subdirPath, 'Decon'), 'dir')
        % Rename the "Decon" folder to match its parent directory name
        movefile(fullfile(subdirPath, 'Decon'), fullfile(subdirPath, [subdirName, '']));
        
        % Move the renamed "Decon" folder to the "Decon" folder in the main folder
        movefile(fullfile(subdirPath, [subdirName, '']), deconFolder);
    end
end
