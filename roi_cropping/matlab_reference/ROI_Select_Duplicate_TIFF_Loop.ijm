for (i = 1; i <= 100; i++) {

    idx = zeroPad2(i);

    open("/Users/theodorecarter/Library/CloudStorage/OneDrive-UniversityofDenver/Asensio CryoData/Theo Carter/EVOS/20260410_Syt9_Insulin_PulseChase/Decon/0min_pma/Deconvoluted/gfp/gfp" + idx + ".tif");

    roiManager("Select", 0);
    run("Duplicate...", "duplicate");
    roiManager("Add");
    run("Make Inverse");
    run("Clear", "stack");

    saveAs("Tiff",
    "/Users/theodorecarter/Library/CloudStorage/OneDrive-UniversityofDenver/Asensio CryoData/Theo Carter/EVOS/20260410_Syt9_Insulin_PulseChase/Decon/0min_pma/Deconvoluted/Cropped/gfp/gfp" + idx + ".tif");

    close();
    close();


    open("/Users/theodorecarter/Library/CloudStorage/OneDrive-UniversityofDenver/Asensio CryoData/Theo Carter/EVOS/20260410_Syt9_Insulin_PulseChase/Decon/0min_pma/Deconvoluted/cy/cy" + idx + ".tif");

    roiManager("Select", 0);
    run("Duplicate...", "duplicate");
    run("Make Inverse");
    run("Clear", "stack");

    saveAs("Tiff",
    "/Users/theodorecarter/Library/CloudStorage/OneDrive-UniversityofDenver/Asensio CryoData/Theo Carter/EVOS/20260410_Syt9_Insulin_PulseChase/Decon/0min_pma/Deconvoluted/Cropped/cy/cy" + idx + ".tif");

    close();
    close();

    roiManager("Delete");
}


// Helper function for leading zeros
function zeroPad2(n) {
    if (n < 10)
        return "0" + n;
    else
        return "" + n;
}