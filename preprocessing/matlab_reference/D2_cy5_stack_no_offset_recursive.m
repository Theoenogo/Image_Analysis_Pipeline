function D2_cy5_stack_no_offset_recursive()

    % Access the global variable mainFolderPath
    global mainFolderPath;

    % Call the function to process "cy" folders in the first subdirectories
    processCYFolders(mainFolderPath, 0);

    function processCYFolders(folderPath, depth)
        % Get a list of all subfolders in the current folder
        subfolderList = dir(folderPath);
        subfolderList = subfolderList([subfolderList.isdir]); % Filter out non-folders
        subfolderList = subfolderList(~ismember({subfolderList.name}, {'.', '..'})); % Exclude '.' and '..'

        % Loop through each subfolder and process the "cy" folders
        for folderIdx = 1:length(subfolderList)
            folderName = subfolderList(folderIdx).name;

            % Check if the folder starts with 'cy'
            if startsWith(folderName, 'cy', 'IgnoreCase', true)
                % Get the current folder path
                cyFolderPath = fullfile(folderPath, folderName);

                % Get a list of all files in the folder with details
                fileListDetails = dir(fullfile(cyFolderPath, '*.TIF'));
                fileListDetails = [fileListDetails; dir(fullfile(cyFolderPath, '*.tiff'))]; % Add lowercase extension

                % Check if any TIFF files are present
                tiffFilesFound = any([fileListDetails.isdir] == 0);

                if ~tiffFilesFound
                    disp(['No TIFF files found in the folder: ' cyFolderPath]);
                else
                    % Initialize an empty cell array to store the stack
                    stack = {};

                    % Loop through each TIFF file and read them to form the stack
                    for i = 1:length(fileListDetails)
                        % Read the current TIFF file
                        currentFile = imread(fullfile(cyFolderPath, fileListDetails(i).name));

                        % Append the current image to the stack
                        stack{i} = currentFile;
                    end

                    % Check if the stack is not empty before saving the stack
                    if ~isempty(stack)
                        % Save the stack as a new multi-page TIFF file
                        stackFileName = [folderName, '.tif'];
                        fullStackFilePath = fullfile(cyFolderPath, stackFileName);

                        % Use imwrite to create a multi-page TIFF file with the stack
                        for i = 1:length(stack)
                            if i == 1
                                imwrite(uint16(stack{i}), fullStackFilePath, 'tif', 'WriteMode', 'overwrite');
                            else
                                imwrite(uint16(stack{i}), fullStackFilePath, 'tif', 'WriteMode', 'append');
                            end
                        end

                        disp(['Stack created for images in folder: ' cyFolderPath]);
                        disp(['Stack saved with folder name as: ' stackFileName]);

                        % Create the "Decon/cy/" path in the parent folder
                        deconCYPath = fullfile(fileparts(cyFolderPath), 'Decon', 'cy');
                        if ~exist(deconCYPath, 'dir')
                            mkdir(deconCYPath);
                        end

                        % Move the newly generated stack to the "Decon/cy/" folder
                        movefile(fullStackFilePath, fullfile(deconCYPath, stackFileName));
                    else
                        disp(['No TIFF files were successfully read in folder: ' cyFolderPath]);
                    end
                end
            elseif depth < 1
                % Recursively process the subfolders if the depth is within the first subdirectories
                processCYFolders(fullfile(folderPath, folderName), depth + 1);
            end
        end
    end
end
