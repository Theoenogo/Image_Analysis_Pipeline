#@ File(label = "Select Background_Subtracted folder", style = "directory") bgSubDir

"""
Combined colocalization analysis for ImageJ (Jython).

Computes per-slice auto-thresholds, Pearson's correlation, and thresholded
Manders' coefficients in a single pass -- no JACoP plugin required.

Key improvements over the separate threshold + JACoP workflow:
  - No GUI windows opened (no show()/close() overhead)
  - No time.sleep() delays
  - No log parsing -- results computed directly
  - Images opened once instead of twice
  - No JACoP plugin dependency

Usage: Run in ImageJ/Fiji via Plugins > Macros > Run...
       Select your Background_Subtracted folder when prompted.

Output:
  coloc_results.csv        -- per-slice Pearson's r, Manders' M1 & M2
  thresholds_multislice.csv -- per-slice thresholds (backward-compatible)
"""

import os
import csv
import math
from ij import IJ, ImagePlus
from ij.process import AutoThresholder

rootDir = bgSubDir.getAbsolutePath()
gfpDir = os.path.join(rootDir, "gfp")
cyDir = os.path.join(rootDir, "cy")


def get_num(filename, prefix):
    """Extract number from filename after prefix."""
    base = os.path.basename(filename)
    s = base[len(prefix):-4]
    try:
        return int(s)
    except ValueError:
        return s


def get_files(folder, prefix):
    """Get sorted list of TIFF files with given prefix."""
    return sorted([f for f in os.listdir(folder)
                   if f.lower().endswith(('.tif', '.tiff'))
                   and f.lower().startswith(prefix.lower())])


def compute_threshold(ip):
    """
    Compute auto-threshold using ImageJ's Default (IsoData) method.

    Uses AutoThresholder directly on the histogram for reliable, correct
    results.  Falls back to mean + 2*stddev if IsoData returns 0.
    """
    histogram = ip.getHistogram()
    thresholder = AutoThresholder()
    thr = thresholder.getThreshold(AutoThresholder.Method.Default, histogram)
    if thr <= 0:
        stats = ip.getStats()
        return stats.mean + 2.0 * stats.stdDev
    return float(thr)


def compute_coloc(ip_a, ip_b, thr_a, thr_b):
    """
    Compute Pearson's correlation and thresholded Manders' coefficients.

    Pearson's r:  standard correlation coefficient over ALL pixels.
    Manders' M1:  sum(A | A>thrA AND B>thrB) / sum(A | A>thrA)
    Manders' M2:  sum(B | A>thrA AND B>thrB) / sum(B | B>thrB)

    Uses FloatProcessor to handle 8-bit, 16-bit, and 32-bit images
    uniformly.
    """
    fp_a = ip_a.convertToFloatProcessor()
    fp_b = ip_b.convertToFloatProcessor()
    pix_a = fp_a.getPixels()
    pix_b = fp_b.getPixels()
    n = len(pix_a)

    # Get means via native Java statistics (fast)
    stats_a = fp_a.getStatistics()
    stats_b = fp_b.getStatistics()
    mean_a = stats_a.mean
    mean_b = stats_b.mean

    # Single pass: Pearson's numerator/denominator + Manders' sums
    num = 0.0
    den_a = 0.0
    den_b = 0.0
    sum_a_above = 0.0
    sum_a_coloc = 0.0
    sum_b_above = 0.0
    sum_b_coloc = 0.0

    for i in range(n):
        a = float(pix_a[i])
        b = float(pix_b[i])

        # Pearson's
        da = a - mean_a
        db = b - mean_b
        num += da * db
        den_a += da * da
        den_b += db * db

        # Thresholded Manders'
        a_above = a > thr_a
        b_above = b > thr_b

        if a_above:
            sum_a_above += a
            if b_above:
                sum_a_coloc += a

        if b_above:
            sum_b_above += b
            if a_above:
                sum_b_coloc += b

    denom = math.sqrt(den_a * den_b)
    pearson = num / denom if denom > 0 else 0.0
    manders_a = sum_a_coloc / sum_a_above if sum_a_above > 0 else 0.0
    manders_b = sum_b_coloc / sum_b_above if sum_b_above > 0 else 0.0

    return pearson, manders_a, manders_b


def main():
    if not os.path.isdir(gfpDir):
        IJ.error("GFP folder not found: " + gfpDir)
        return
    if not os.path.isdir(cyDir):
        IJ.error("CY folder not found: " + cyDir)
        return

    results_csv = os.path.join(rootDir, "coloc_results.csv")
    thresholds_csv = os.path.join(rootDir, "thresholds_multislice.csv")

    print("GFP folder: " + gfpDir)
    print("CY folder:  " + cyDir)
    print("Results:     " + results_csv)

    gfp_files = get_files(gfpDir, "gfp")
    cy_files = get_files(cyDir, "cy")
    print("Found {} GFP and {} CY files".format(len(gfp_files), len(cy_files)))

    # Write CSV headers
    with open(results_csv, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(["gfp_image", "cy_image", "slice",
                          "gfp_threshold", "cy_threshold",
                          "pearson", "mandersA", "mandersB"])

    with open(thresholds_csv, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(["gfp_image", "cy_image", "slice",
                          "gfp_threshold", "cy_threshold"])

    total_slices = 0
    total_pairs = 0

    for gfp_file in gfp_files:
        gfp_num = get_num(gfp_file, "gfp")
        cy_file = next((f for f in cy_files if get_num(f, "cy") == gfp_num), None)

        if not cy_file:
            print("No matching CY file for " + gfp_file)
            continue

        print("\nProcessing: {} + {}".format(gfp_file, cy_file))

        gfp_stack = IJ.openImage(os.path.join(gfpDir, gfp_file))
        cy_stack = IJ.openImage(os.path.join(cyDir, cy_file))

        if gfp_stack is None or cy_stack is None:
            print("Error: could not open one or both images")
            if gfp_stack:
                gfp_stack.close()
            if cy_stack:
                cy_stack.close()
            continue

        nslices = min(gfp_stack.getNSlices(), cy_stack.getNSlices())
        print("  {} slices".format(nslices))

        result_rows = []
        threshold_rows = []

        for s in range(1, nslices + 1):
            ip_g = gfp_stack.getStack().getProcessor(s)
            ip_c = cy_stack.getStack().getProcessor(s)

            thr_g = compute_threshold(ip_g)
            thr_c = compute_threshold(ip_c)

            pearson, m_a, m_b = compute_coloc(ip_g, ip_c, thr_g, thr_c)

            print("  Slice {:3d}: thr_g={:7.1f}  thr_c={:7.1f}  r={:.4f}  M1={:.4f}  M2={:.4f}".format(
                s, thr_g, thr_c, pearson, m_a, m_b))

            result_rows.append([gfp_file, cy_file, s,
                                "{:.2f}".format(thr_g), "{:.2f}".format(thr_c),
                                "{:.6f}".format(pearson),
                                "{:.6f}".format(m_a),
                                "{:.6f}".format(m_b)])
            threshold_rows.append([gfp_file, cy_file, s,
                                   "{:.2f}".format(thr_g), "{:.2f}".format(thr_c)])

        # Batch-write all slices for this pair
        with open(results_csv, 'a') as f:
            csv.writer(f).writerows(result_rows)
        with open(thresholds_csv, 'a') as f:
            csv.writer(f).writerows(threshold_rows)

        total_slices += nslices
        total_pairs += 1

        gfp_stack.close()
        cy_stack.close()

    print("\n========================================")
    print("Complete!")
    print("Processed {} image pairs, {} total slices".format(total_pairs, total_slices))
    print("Results:    " + results_csv)
    print("Thresholds: " + thresholds_csv)


main()
