// =====================================
// Select folders and ROI zip
// =====================================

showMessage("Select GFP Folder");
gfpDir = getDirectory("Choose GFP folder");
if (gfpDir == "") exit("Cancelled.");

showMessage("Select CY Folder");
cyDir = getDirectory("Choose CY folder");
if (cyDir == "") exit("Cancelled.");

roiZipPath = File.openDialog("Choose the single roi.zip file");
if (roiZipPath == "") exit("Cancelled.");


// =====================================
// Prompt subtraction multipliers
// =====================================

gfpMultiplier = getNumber("Enter GFP subtraction multiplier:", 1.25);
cyMultiplier  = getNumber("Enter CY subtraction multiplier:", 1.25);


// =====================================
// Output folders
// =====================================

parentFolder = getSharedParent(gfpDir, cyDir);
if (parentFolder == "")
    exit("GFP and CY folders are not in the same parent directory.");

outputRoot = parentFolder + "Background_Subtracted" + File.separator;
gfpOutput  = outputRoot + "gfp" + File.separator;
cyOutput   = outputRoot + "cy" + File.separator;

File.makeDirectory(outputRoot);
File.makeDirectory(gfpOutput);
File.makeDirectory(cyOutput);


// =====================================
// Build sorted image lists
// =====================================

gfpFiles = getSortedTiffList(gfpDir);
cyFiles  = getSortedTiffList(cyDir);

if (gfpFiles.length == 0)
    exit("No GFP TIFF files found.");

if (cyFiles.length == 0)
    exit("No CY TIFF files found.");

if (gfpFiles.length != cyFiles.length)
    exit("GFP and CY folders do not contain the same number of TIFF files.");


// =====================================
// Load ROI zip
// =====================================

run("ROI Manager...");
roiManager("Reset");
roiManager("Open", roiZipPath);

roiCount = roiManager("Count");

if (roiCount == 0)
    exit("No ROIs found in selected roi.zip.");

if (roiCount < gfpFiles.length)
    exit("Not enough ROIs in roi.zip for the number of images.\nROIs: " + roiCount + "\nImages: " + gfpFiles.length);

print("Loaded ROI count: " + roiCount);
print("GFP image count: " + gfpFiles.length);
print("CY image count: " + cyFiles.length);


// =====================================
// Main processing loop
// ROI 0 -> image 1
// ROI 1 -> image 2
// etc.
// =====================================

for (i = 0; i < gfpFiles.length; i++) {

    gfpName = gfpFiles[i];
    cyName  = cyFiles[i];

    print("Processing image pair:");
    print("  GFP: " + gfpName);
    print("  CY:  " + cyName);
    print("  ROI index: " + i);

    processImageOriginalLogic(gfpDir + gfpName, gfpOutput, gfpMultiplier, "GFP", i);
    processImageOriginalLogic(cyDir + cyName,  cyOutput,  cyMultiplier,  "CY",  i);
}

showMessage("Processing complete.");


// =====================================
// Process one image using original logic
// =====================================

function processImageOriginalLogic(path, outputDir, multiplier, channelLabel, roiIndex) {

    fileName = File.getName(path);

    open(path);

    if (nSlices < 1) {
        print(channelLabel + " | Skipping file due to loading error: " + fileName);
        if (nImages > 0) close();
        return;
    }

    // Apply ROI to the opened image
    roiManager("Select", roiIndex);
    print(channelLabel + " | Selecting ROI index: " + roiIndex);

    run("Clear Results");
    run("Measure");

    meanIntensity = getResult("Mean", 0);
    print(channelLabel + " | Measured mean intensity: " + meanIntensity);

    // ORIGINAL LOGIC: one subtraction value per image
    subtractValue = meanIntensity * multiplier;
    if (subtractValue < 100) subtractValue = 100;
    if (subtractValue > 5000) subtractValue = 5000;

    print(channelLabel + " | Calculated subtraction value: " + subtractValue);

    if (nSlices > 1) {
        for (s = 1; s <= nSlices; s++) {
            setSlice(s);
            print(channelLabel + " | Subtracting value: " + subtractValue + " from slice: " + s);
            run("Subtract...", "value=" + subtractValue);
        }
    } else {
        print(channelLabel + " | Subtracting value: " + subtractValue + " from single slice");
        run("Subtract...", "value=" + subtractValue);
    }

    saveAs("Tiff", outputDir + fileName);
    close();
}


// =====================================
// Sort TIFF list numerically
// =====================================

function getSortedTiffList(folder) {

    raw = getFileList(folder);

    count = 0;
    for (i = 0; i < raw.length; i++) {
        if (endsWith(raw[i], ".tif") || endsWith(raw[i], ".tiff"))
            count++;
    }

    list = newArray(count);
    nums = newArray(count);

    j = 0;
    for (i = 0; i < raw.length; i++) {
        fileName = raw[i];

        if (!endsWith(fileName, ".tif") && !endsWith(fileName, ".tiff"))
            continue;

        list[j] = fileName;
        nums[j] = getNumericIndex(fileName);
        j++;
    }

    for (a = 0; a < count - 1; a++) {
        for (b = a + 1; b < count; b++) {
            if (nums[a] > nums[b]) {
                tempNum = nums[a];
                nums[a] = nums[b];
                nums[b] = tempNum;

                tempName = list[a];
                list[a] = list[b];
                list[b] = tempName;
            }
        }
    }

    return list;
}


// =====================================
// Extract first numeric block from filename
// =====================================

function getNumericIndex(name) {

    s = replace(name, ".tif", "");
    s = replace(s, ".tiff", "");

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

    if (numStr == "")
        return 999999;

    return 0 + numStr;
}


// =====================================
// Shared parent helper
// =====================================

function getSharedParent(path1, path2) {

    p1 = stripTrailingSeparator(path1);
    p2 = stripTrailingSeparator(path2);

    parent1 = substring(p1, 0, lastIndexOf(p1, File.separator) + 1);
    parent2 = substring(p2, 0, lastIndexOf(p2, File.separator) + 1);

    if (parent1 == parent2)
        return parent1;

    return "";
}


// =====================================
// Strip trailing separator
// =====================================

function stripTrailingSeparator(path) {
    if (endsWith(path, File.separator))
        return substring(path, 0, lengthOf(path) - 1);
    return path;
}