%% Step 1: Prompt user for experiment base directory

baseFolderPath = uigetdir('', 'Select the experiment base folder (contains gfp/cy folders)');

if isequal(baseFolderPath,0)
    error('No folder selected. Script cancelled.');
end

%% Step 2: Prompt for channel folder names

prompt = {'Enter image channel folder names (comma separated):'};
dlgtitle = 'Channel Folder Names';
dims = [1 50];
defaultInput = {'gfp,cy'};
answer = inputdlg(prompt,dlgtitle,dims,defaultInput);

if isempty(answer)
    error('Operation cancelled.');
end

imageFolderNames = strtrim(strsplit(answer{1},','));

%% Step 3: Prompt for ROI folder name

roiAnswer = inputdlg('Enter ROI folder name:','ROI Folder',[1 50],{'roi'});

if isempty(roiAnswer)
    error('Operation cancelled.');
end

roiFolderName = roiAnswer{1};

%% Confirm operation

choice = questdlg(sprintf('Process folders inside:\n%s\n\nContinue?',baseFolderPath), ...
    'Confirm Operation','Yes','Cancel','Cancel');

if ~strcmp(choice,'Yes')
    error('Operation cancelled.');
end


%% Step 4: Find all folders recursively

allFolders = dir(fullfile(baseFolderPath,'**'));
allFolders = allFolders([allFolders.isdir]);


%% Step 5: Process directories

for f = 1:numel(allFolders)

    currentFolder = allFolders(f).folder;
    currentFolderName = allFolders(f).name;

    if strcmp(currentFolderName,'.') || strcmp(currentFolderName,'..')
        continue;
    end


    %% Check if folder matches channel folders

    if any(strcmp(currentFolderName,imageFolderNames))

        imageFolderPath = fullfile(currentFolder,currentFolderName);

        parentFolderPath = fileparts(imageFolderPath);

        roiZipFolderPath = fullfile(parentFolderPath,roiFolderName);


        %% Verify ROI folder exists

        if ~exist(roiZipFolderPath,'dir')

            fprintf('Skipping %s (ROI folder missing)\n',imageFolderPath);

            continue;

        end


        %% Get image + ROI zip lists

        imageFiles = dir(fullfile(imageFolderPath,'*.tif'));

        zipFiles = dir(fullfile(roiZipFolderPath,'*.zip'));


        %% Process ROI zip files

        for i = 1:numel(zipFiles)

            tempFolder = fullfile(parentFolderPath,'temp_roi_extract');

            unzip(fullfile(roiZipFolderPath,zipFiles(i).name),tempFolder);

            roiFiles = dir(fullfile(tempFolder,'*.roi'));

            numRoiFiles = numel(roiFiles);


            %% Duplicate images if multiple ROIs

            if i <= numel(imageFiles) && numRoiFiles >= 2

                imageFileName = imageFiles(i).name;

                imageFilePath = fullfile(imageFolderPath,imageFileName);

                for j = 2:numRoiFiles

                    newImageFileName = sprintf('%s_%d.tif',imageFileName(1:end-4),j);

                    copyfile(imageFilePath,fullfile(imageFolderPath,newImageFileName));

                end

            elseif numRoiFiles < 2

                fprintf('Skipping duplication (single ROI only)\n');

            else

                fprintf('Mismatch between ROI and image count\n');

                break;

            end


            %% Remove temp extraction folder

            rmdir(tempFolder,'s');

        end


        %% Backup originals

        backupFolderPath = fullfile(imageFolderPath,'backup');

        if ~exist(backupFolderPath,'dir')

            mkdir(backupFolderPath);

        end


        for i = 1:numel(imageFiles)

            copyfile(fullfile(imageFolderPath,imageFiles(i).name),backupFolderPath);

        end


        %% Rename images sequentially

        renamedFiles = dir(fullfile(imageFolderPath,'*.tif'));

        for i = 1:numel(renamedFiles)

            newImageFileName = sprintf('%s%d.tif',currentFolderName,i);

            movefile(fullfile(imageFolderPath,renamedFiles(i).name), ...
                     fullfile(imageFolderPath,newImageFileName));

        end


        %% Create Cropped folder structure

        croppedFolderPath = fullfile(parentFolderPath,'Cropped');

        if ~exist(croppedFolderPath,'dir')

            mkdir(croppedFolderPath);

        end


        updatedRenamedFiles = dir(fullfile(imageFolderPath,'*.tif'));

        for i = 1:numel(updatedRenamedFiles)

            [~,imageName,~] = fileparts(updatedRenamedFiles(i).name);

            mkdir(fullfile(croppedFolderPath,imageName));

        end


        fprintf('Finished processing folder: %s\n',imageFolderPath);

    end

end


fprintf('\nProcessing complete.\n');