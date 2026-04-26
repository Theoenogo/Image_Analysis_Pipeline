function B_deleteScanProtocolFiles()

    % Access the global variable mainFolderPath
    global mainFolderPath;

    % Call the function with the given folder path
    deleteScanProtocolFilesRecursive(mainFolderPath);

    function deleteScanProtocolFilesRecursive(folderPath)
        files = dir(folderPath);
        for i = 1:numel(files)
            if ~strcmp(files(i).name, '.') && ~strcmp(files(i).name, '..')
                filePath = fullfile(folderPath, files(i).name);

                if files(i).isdir
                    % Recursively search subdirectories
                    deleteScanProtocolFilesRecursive(filePath);
                else
                    [~, ~, fileExt] = fileparts(filePath);
                    if strcmpi(fileExt, '.scanprotocol')
                        delete(filePath);
                        disp(['Deleted: ', filePath]);
                    end
                end
            end
        end
    end

end
