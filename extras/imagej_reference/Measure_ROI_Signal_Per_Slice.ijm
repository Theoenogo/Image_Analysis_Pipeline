// Prompt user to select input folder for cy .tif files
cyInputDir = getDirectory("Choose a folder containing the cy .tif files");
if (cyInputDir == "") {
    showMessage("No folder selected. Process cancelled.");
    exit();
}

// Determine parent directory of cy folder
parentDir = substring(cyInputDir, 0, lastIndexOf(cyInputDir, File.separator));

// Create Results folder in parent directory
outputDir = parentDir + "Results" + File.separator;
File.makeDirectory(outputDir);

print("Saving CSV results to: " + outputDir);

// Get full file list
fileList = getFileList(cyInputDir);

// Build filtered TIFF list (cy*.tif only)
cyFiles = newArray();
indexList = newArray();

for (i = 0; i < fileList.length; i++) {

    name = fileList[i];

    if (endsWith(name.toLowerCase(), ".tif") &&
        startsWith(name.toLowerCase(), "cy")) {

        cyFiles = Array.concat(cyFiles, name);
        indexList = Array.concat(indexList, getNumericIndex(name));
    }
}

// Sort numerically (cy01 → cy02 → cy03)
for (a = 0; a < cyFiles.length - 1; a++) {

    for (b = a + 1; b < cyFiles.length; b++) {

        if (indexList[a] > indexList[b]) {

            tempNum = indexList[a];
            indexList[a] = indexList[b];
            indexList[b] = tempNum;

            tempName = cyFiles[a];
            cyFiles[a] = cyFiles[b];
            cyFiles[b] = tempName;
        }
    }
}

print("Number of cy TIFF files detected: " + cyFiles.length);

// Confirm ROI Manager loaded
numROIs = roiManager("count");

if (numROIs == 0) {
    showMessage("No ROIs loaded in the ROI Manager.");
    exit();
}

// Process images
for (i = 0; i < cyFiles.length; i++) {

    cyFileName = cyFiles[i];

    print("Processing: " + cyFileName);

    open(cyInputDir + cyFileName);

    if (isOpen(cyFileName)) {

        stackSize = nSlices;

        roiIndex = i % numROIs;
        roiManager("Select", roiIndex);

        print("Using ROI index: " + roiIndex);

        run("Clear Results");

        for (slice = 1; slice <= stackSize; slice++) {

            setSlice(slice);

            run("Duplicate...", "title=TempSlice");

            selectWindow("TempSlice");

            roiManager("Select", roiIndex);

            run("Measure");

            close();
        }

        csvFileName = replace(cyFileName, ".tif", ".csv");

        saveAs("Results", outputDir + csvFileName);

        print("Saved CSV: " + csvFileName);

        run("Clear Results");

        close();

    } else {

        print("Failed to open: " + cyFileName);

    }
}

showMessage("Processing completed.\nResults saved to:\n" + outputDir);


// Helper: extract numeric portion from filename
function getNumericIndex(name) {

    s = replace(name, ".tif", "");

    numStr = "";

    foundDigit = false;

    for (k = 0; k < lengthOf(s); k++) {

        ch = substring(s, k, k + 1);

        if (ch >= "0" && ch <= "9") {

            numStr = numStr + ch;
            foundDigit = true;

        } else if (foundDigit) {

            break;

        }
    }

    if (numStr == "") return 999999;

    return 0 + numStr;
}
