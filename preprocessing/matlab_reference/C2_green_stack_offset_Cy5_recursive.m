function C2_green_stack_offset_recursive()

    % Access the global variable mainFolderPath
    global mainFolderPath;

    % Call the function to process "gfp" folders in the first subdirectories
    processGFPFolders(mainFolderPath, 0);

    function processGFPFolders(folderPath, depth)
        % Get a list of all subfolders in the current folder
        subfolderList = dir(folderPath);
        subfolderList = subfolderList([subfolderList.isdir]); % Filter out non-folders
        subfolderList = subfolderList(~ismember({subfolderList.name}, {'.', '..'})); % Exclude '.' and '..'

        % Loop through each subfolder and process the "gfp" folders
        for folderIdx = 1:length(subfolderList)
            folderName = subfolderList(folderIdx).name;

            % Check if the folder starts with 'gfp'
            if startsWith(folderName, 'gfp', 'IgnoreCase', true)
                % Get the current folder path
                gfpFolderPath = fullfile(folderPath, folderName);

                % Get a list of all files in the folder with details
                fileListDetails = dir(fullfile(gfpFolderPath, '*.TIF'));
                fileListDetails = [fileListDetails; dir(fullfile(gfpFolderPath, '*.tiff'))]; % Add lowercase extension

                % Check if any TIFF files are present
                tiffFilesFound = any([fileListDetails.isdir] == 0);

                if ~tiffFilesFound
                    disp(['No TIFF files found in the folder: ' gfpFolderPath]);
                else
                    % Initialize an empty cell array to store the stack
                    stack = {};

                    % Loop through each TIFF file and read them to form the stack
                    for i = 1:length(fileListDetails)
                        % Read the current TIFF file
                        currentFile = imread(fullfile(gfpFolderPath, fileListDetails(i).name));

                        % Apply XY offset to the current image using affine transformation
                        xOffset = 5;
                        yOffset = -2;
                        tform = affine2d([1 0 0; 0 1 0; xOffset yOffset 1]);
                        offsetImage = imwarp(currentFile, tform, 'OutputView', imref2d(size(currentFile)));

                        % Append the current offset image to the stack
                        stack{i} = offsetImage;
                    end

                    % Check if the stack is not empty before saving the offset image
                    if ~isempty(stack)
                        % Save the offset image as a new multi-page TIFF file
                        offsetFileName = [folderName, '.tif'];
                        fullOffsetFilePath = fullfile(gfpFolderPath, offsetFileName);

                        % Use imwrite to create a multi-page TIFF file with the offset stack
                        for i = 1:length(stack)
                            if i == 1
                                imwrite(uint16(stack{i}), fullOffsetFilePath, 'tif', 'WriteMode', 'overwrite');
                            else
                                imwrite(uint16(stack{i}), fullOffsetFilePath, 'tif', 'WriteMode', 'append');
                            end
                        end

                        disp(['XY offset applied to the images in folder: ' gfpFolderPath]);
                        disp(['Offset image saved with folder name as: ' offsetFileName]);

                        % Create the "Decon/gfp/" path in the parent folder
                        deconGFPPath = fullfile(fileparts(gfpFolderPath), 'Decon', 'gfp');
                        if ~exist(deconGFPPath, 'dir')
                            mkdir(deconGFPPath);
                        end

                        % Move the newly generated offset image to the "Decon/gfp/" folder
                        movefile(fullOffsetFilePath, fullfile(deconGFPPath, offsetFileName));
                    else
                        disp(['No TIFF files were successfully read in folder: ' gfpFolderPath]);
                    end
                end
            elseif depth < 1
                % Recursively process the subfolders if the depth is within the first subdirectories
                processGFPFolders(fullfile(folderPath, folderName), depth + 1);
            end
        end
    end
end
