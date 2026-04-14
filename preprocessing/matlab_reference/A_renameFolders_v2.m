function A_renameFolders_v2()

    % Access the global variable mainFolderPath
    global mainFolderPath;
    
    function renameFoldersInDirectories(pattern, directoryPath)
        % Get a list of all subdirectories in the specified directory
        directories = dir(directoryPath);
        directories = directories([directories.isdir]);

        % Loop through each subdirectory and rename folders if they match the pattern
        for i = 1:numel(directories)
            directoryName = directories(i).name;
            if ~strcmp(directoryName, '.') && ~strcmp(directoryName, '..')
                subDirPath = fullfile(directoryPath, directoryName);

                % Get a list of all folders in the current subdirectory
                folders = dir(subDirPath);
                folders = folders([folders.isdir]);

                % Loop through each folder and rename if it matches the pattern
                for j = 1:numel(folders)
                    folderName = folders(j).name;
                    if ~strcmp(folderName, '.') && ~strcmp(folderName, '..') && ~isempty(regexp(folderName, pattern, 'once'))
                        % Extract the part of the folder name before the dot
                        newName = strtok(folderName, '.');

                        % Rename the folder
                        movefile(fullfile(subDirPath, folderName), fullfile(subDirPath, newName));
                    end
                end
            end
        end
    end

    % Define the regular expression patterns for matching folders
    gfpPattern = 'gfp\d+\.[\dA-Za-z]+';
    cyPattern = 'cy\d+\.[\dA-Za-z]+';
    rfpPattern = 'rfp\d+\.[\dA-Za-z]+';

    % Call the function to rename folders in all subdirectories of the specified directory path
renameFoldersInDirectories(gfpPattern, mainFolderPath);
renameFoldersInDirectories(cyPattern, mainFolderPath);
renameFoldersInDirectories(rfpPattern, mainFolderPath);

end
