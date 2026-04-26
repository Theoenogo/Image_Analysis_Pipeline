// =====================================
// Prompt user before folder selection
// =====================================

showMessage("Select GFP Folder",
    "In the next dialog, choose the folder containing your GFP images.");

pathGFP = getDirectory("Choose GFP image folder");

if (pathGFP == "")
    exit("GFP folder selection cancelled.");

showMessage("Select CY Folder",
    "In the next dialog, choose the folder containing your CY images.");

pathCY = getDirectory("Choose CY image folder");

if (pathCY == "")
    exit("CY folder selection cancelled.");

showMessage("Select ROI Source Folder",
    "In the next dialog, choose the folder containing your existing ROI zip files (e.g. roi_original).");

pathROISource = getDirectory("Choose ROI source folder");

if (pathROISource == "")
    exit("ROI source folder selection cancelled.");


// =====================================
// Determine parent folder and ROI folder
// =====================================

parentGFP = getParentFolder(pathGFP);
parentCY  = getParentFolder(pathCY);

if (parentGFP != parentCY) {
    exit("The GFP and CY folders must be inside the same parent directory.");
}

roiFolder = parentGFP + "roi" + File.separator;

if (!File.exists(roiFolder)) {
    File.makeDirectory(roiFolder);
}

print("GFP folder: " + pathGFP);
print("CY folder: " + pathCY);
print("ROI source folder: " + pathROISource);
print("ROI save folder: " + roiFolder);


// =====================================
// Get and sort GFP file list naturally
// =====================================

rawList = getFileList(pathGFP);

// Keep only gfp tif/tiff files
count = 0;
for (i = 0; i < rawList.length; i++) {
    f = rawList[i];
    if ((!endsWith(f, ".tif")) && (!endsWith(f, ".tiff")))
        continue;
    if (!startsWith(f, "gfp"))
        continue;
    count++;
}

if (count == 0)
    exit("No GFP TIFF files found.");

fileList = newArray(count);
fileNums = newArray(count);

j = 0;
for (i = 0; i < rawList.length; i++) {
    f = rawList[i];
    if ((!endsWith(f, ".tif")) && (!endsWith(f, ".tiff")))
        continue;
    if (!startsWith(f, "gfp"))
        continue;

    fileList[j] = f;
    fileNums[j] = getNumericIndex(f);
    j++;
}

// Sort by numeric index
for (a = 0; a < count - 1; a++) {
    for (b = a + 1; b < count; b++) {
        if (fileNums[a] > fileNums[b]) {

            tempNum = fileNums[a];
            fileNums[a] = fileNums[b];
            fileNums[b] = tempNum;

            tempName = fileList[a];
            fileList[a] = fileList[b];
            fileList[b] = tempName;
        }
    }
}

print("Found " + fileList.length + " GFP images after filtering/sorting.");


// =====================================
// Open ROI Manager once
// =====================================

run("ROI Manager...");


// =====================================
// Loop through GFP images
// =====================================

for (i = 0; i < fileList.length; i++) {

    filename = fileList[i];
    imageIndex = fileNums[i];

    if (imageIndex < 0) {
        print("Skipping file with no numeric index: " + filename);
        continue;
    }

    cyName = replace(filename, "gfp", "cy");

    if (!File.exists(pathCY + cyName)) {
        print("Skipping missing CY pair for: " + filename);
        continue;
    }

    // =====================================
    // Create padded ROI filename (01.zip, 02.zip, ...)
    // =====================================

    if (imageIndex < 10)
        roiName = "0" + imageIndex + ".zip";
    else
        roiName = "" + imageIndex + ".zip";

    roiPath = roiFolder + roiName;

    // Reset ROI Manager
    roiManager("Reset");

    // =====================================
    // Open channels and merge
    // =====================================

    open(pathCY + cyName);
    open(pathGFP + filename);

    run("Merge Channels...",
        "c1=[" + cyName + "] c2=[" + filename + "] create");

    mergedTitle = "merged_" + filename;

    selectWindow("Composite");
    rename(mergedTitle);

    // =====================================
    // If an ROI zip exists in the source
    // folder, load it into the manager
    // =====================================

    roiSourcePath = pathROISource + roiName;

    if (File.exists(roiSourcePath)) {
        roiManager("Open", roiSourcePath);
        print("Loaded existing ROIs from: " + roiSourcePath);

        waitForUser(
            "Editing ROIs for:\n" + filename +
            "\n\nExisting ROIs have been loaded. Add or remove as needed." +
            "\n\nClick OK when finished.\nROIs will be renamed and saved automatically."
        );

    } else {

        // =====================================
        // No existing ROIs — draw from scratch
        // =====================================

        waitForUser(
            "Draw ROIs for:\n" + filename +
            "\n\nClick OK when finished.\nROIs will be renamed and saved automatically."
        );
    }

    // =====================================
    // Rename ROIs
    // Example:
    // 1-1, 1-2
    // 10-1, 10-2
    // =====================================

    roiTotal = roiManager("count");

    if (roiTotal > 0) {

        for (r = 0; r < roiTotal; r++) {

            roiManager("Select", r);

            roiLabel = "" + imageIndex + "-" + (r + 1);

            roiManager("Rename", roiLabel);
        }

        roiManager("Save", roiPath);

        print("Saved ROI set: " + roiPath);

    } else {

        print("No ROIs drawn for: " + filename);
    }

    // =====================================
    // Close windows safely
    // =====================================

    if (isOpen(cyName)) {
        selectWindow(cyName);
        close();
    }

    if (isOpen(filename)) {
        selectWindow(filename);
        close();
    }

    if (isOpen(mergedTitle)) {
        selectWindow(mergedTitle);
        close();
    }
}

print("Finished processing all images.");


// =====================================
// Helper: parent folder
// =====================================

function getParentFolder(path) {

    p = path;

    if (endsWith(p, File.separator))
        p = substring(p, 0, lengthOf(p) - 1);

    lastSep = lastIndexOf(p, File.separator);

    if (lastSep == -1)
        return "";

    return substring(p, 0, lastSep + 1);
}


// =====================================
// Helper: numeric index extractor
// Example:
// gfp1_decon.tif   -> 1
// gfp10.tif        -> 10
// gfp03_something  -> 3
// =====================================

function getNumericIndex(name) {

    s = replace(name, "gfp", "");
    s = replace(s, ".tif", "");
    s = replace(s, ".tiff", "");

    numStr = "";

    for (k = 0; k < lengthOf(s); k++) {
        ch = substring(s, k, k + 1);

        if (ch >= "0" && ch <= "9")
            numStr = numStr + ch;
        else
            break;
    }

    if (numStr == "")
        return -1;

    return 0 + numStr;
}
