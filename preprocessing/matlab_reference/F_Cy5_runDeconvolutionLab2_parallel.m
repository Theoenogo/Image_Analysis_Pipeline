function runDeconvolutionJob(psfImagePath, baseFolderPath)

    % Load the PSF
    infoPSF = imfinfo(psfImagePath);
    numSlicesPSF = numel(infoPSF);
    psf = zeros(infoPSF(1).Height, infoPSF(1).Width, numSlicesPSF);

    for k = 1:numSlicesPSF
        psf(:,:,k) = imread(psfImagePath, 'Index', k);
    end

    % Convert the PSF to type double for mathematical operations
    psf = double(psf);

    % Find all 'cy' folders recursively within the base folder
    cyDirs = dir(fullfile(baseFolderPath, '**', 'cy'));

    for dirIdx = 1:length(cyDirs)
        currentCyDir = fullfile(cyDirs(dirIdx).folder, cyDirs(dirIdx).name);

        % Process all images within this 'cy' directory
        processCyDirectory(currentCyDir, psf);

        % Move deconvoluted files
        moveDeconFiles(currentCyDir);
    end
    
     % Close MATLAB instance after the script is done
    exit;
end


function processCyDirectory(folderPath, psf)
    % List all the images in the folder
    imageFiles = dir(fullfile(folderPath, '*.tif'));

    for i = 1:length(imageFiles)
        inputImagePath = fullfile(folderPath, imageFiles(i).name);

        % Load the entire image stack
        info = imfinfo(inputImagePath);
        numSlices = numel(info);
        inputImage = zeros(info(1).Height, info(1).Width, numSlices);

        for k = 1:numSlices
            inputImage(:,:,k) = imread(inputImagePath, 'Index', k);
        end

        % Convert images to type double for mathematical operations
        inputImage = double(inputImage);

        % Run the Job function to perform deconvolution
        result = Job(inputImage, psf);

        % Check the minimum and maximum of the result
        minVal = min(result(:));
        maxVal = max(result(:));

        % Rescale if needed
        if minVal < 0 || maxVal > 65535
            result = (result - minVal) / (maxVal - minVal) * 65535;
        end

        % Convert result to uint16
        result16bit = uint16(result);

        % Save the result
        outputPath = fullfile(folderPath, [imageFiles(i).name(1:end-4) '_decon.tif']);
        for k = 1:size(result16bit, 3)
            if k == 1
                imwrite(result16bit(:,:,k), outputPath, 'WriteMode', 'overwrite', 'Compression', 'none');
            else
                imwrite(result16bit(:,:,k), outputPath, 'WriteMode', 'append', 'Compression', 'none');
            end
        end

        % Close Java windows
        closeImageJWindows();
        closeSpecificJavaWindow('Monitor of Matlab RL');

        % Display a message
        disp(['Deconvolved image saved to: ' outputPath]);
    end
end

function result = Job(image, psf)
    % Add DeconvolutionLab_2.jar to the Java path
    javaaddpath([matlabroot filesep 'java' filesep 'DeconvolutionLab_2.jar']);
    
    % Run the deconvolution
    result = DL2.RL(image, psf, 30.0000 , '-out stack short');
end

function moveDeconFiles(currentCyDir)
    % Define the path for the new directory structure
    parentDir = fileparts(currentCyDir);
    grandParentDir = fileparts(parentDir);
    newDeconDir = fullfile(grandParentDir, 'Deconvoluted', 'cy');

    if ~exist(newDeconDir, 'dir')
        mkdir(newDeconDir);
    end

    % Move all _decon files
    deconFiles = dir(fullfile(currentCyDir, '*_decon.tif'));
    for i = 1:length(deconFiles)
        movefile(fullfile(deconFiles(i).folder, deconFiles(i).name), newDeconDir);
    end

    disp(['Moved deconvoluted images to: ' newDeconDir]);
end


function closeImageJWindows()
    % Close ImageJ windows
    ij.WindowManager.closeAllWindows();
end

function closeSpecificJavaWindow(windowTitle)
    % Close a specific Java window based on its title
    frames = java.awt.Frame.getFrames();

    for idx = 1:numel(frames)
        frame = frames(idx);
        if strcmp(frame.getTitle(), windowTitle)
            frame.dispose();
        end
    end
end
