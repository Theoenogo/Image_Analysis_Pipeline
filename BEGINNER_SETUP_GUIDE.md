# Beginner Setup Guide — Image Analysis Pipeline

> **This guide is for you if:** you have never used a terminal, command prompt, or written any code before. Every step is explained in plain language. Follow it top to bottom and you will have the pipeline running.

---

## Table of Contents

1. [What does this pipeline do?](#1-what-does-this-pipeline-do)
2. [What you need before you start](#2-what-you-need-before-you-start)
3. [How to open a terminal](#3-how-to-open-a-terminal)
4. [Install Python 3.12](#4-install-python-312)
5. [Install Git](#5-install-git)
6. [Download the pipeline code](#6-download-the-pipeline-code)
7. [Create a virtual environment](#7-create-a-virtual-environment)
8. [Install the pipeline's dependencies](#8-install-the-pipelines-dependencies)
9. [Running the pipeline — step by step](#9-running-the-pipeline--step-by-step)
10. [Every time you come back: the daily checklist](#10-every-time-you-come-back-the-daily-checklist)
11. [Troubleshooting common problems](#11-troubleshooting-common-problems)

---

## 1. What does this pipeline do?

This is a set of Python programs that process microscopy images of cells. You give it raw image files from the microscope and it automatically:

1. **Cleans and stacks** the raw images
2. **Detects and draws outlines** (called ROIs) around individual cells
3. *(You manually check and fix the outlines in ImageJ — the only manual step)*
4. **Crops** each cell out into its own image file
5. **Removes background light** from each image
6. **Measures how much two fluorescent signals overlap** inside each cell (colocalization)

Think of it like an assembly line: raw images go in one end, measurements come out the other.

---

## 2. What you need before you start

You need to install two free programs before anything else:

| Program | What it is | Why you need it |
|---------|-----------|-----------------|
| **Python 3.12** | The language the pipeline is written in | Runs the scripts |
| **Git** | A program for downloading code | Downloads the pipeline to your computer |

Both are free. Both installation steps are covered below.

---

## 3. How to open a terminal

A terminal (also called "command prompt" or "command line") is a text-based window where you type instructions directly to your computer. It looks plain and old-fashioned but it is very powerful. You do not need to understand how it works — just follow the commands exactly as written.

### On a Mac

1. Press **Command (⌘) + Space** at the same time. A search bar called Spotlight appears.
2. Type `Terminal` and press **Enter**.
3. A window opens with some text ending in a `%` or `$` symbol. That symbol is the "prompt" — it just means the computer is ready for you to type. **You do not type the `$` or `%` yourself.**

Alternatively: open **Finder** → **Applications** → **Utilities** → double-click **Terminal**.

### On Windows

1. Click the **Start button** (the Windows logo in the bottom-left corner).
2. Type `PowerShell` and press **Enter**.
3. A blue or dark window opens with text ending in `>`. That `>` is the prompt — **do not type it yourself.**

> **Note — Windows users:** this guide uses **PowerShell** throughout. If you open the wrong one by accident, close it and try again.

### Reading commands in this guide

Commands in this guide look like this:

```
python --version
```

That means: click inside the terminal window, type exactly what is shown, then press **Enter**. The computer will respond with some text below your line.

---

## 4. Install Python 3.12

Python is the programming language the pipeline is written in. You need version **3.12 specifically** — other versions may not work.

### Check if Python 3.12 is already installed

Open a terminal and type:

```
python3.12 --version
```

- If you see `Python 3.12.x` — **you already have it, skip to Section 5.**
- If you see an error or a different version number — follow the steps below.

### Install Python 3.12 on a Mac

1. Open a web browser and go to: **https://www.python.org/downloads/**
2. Look for the download button for **Python 3.12** (you may need to scroll down past the latest version to find 3.12 specifically — look for a section called "Looking for a specific release?" and find any `3.12.x` version).
3. Click the link for `macOS 64-bit universal2 installer`.
4. A file called something like `python-3.12.x-macos11.pkg` downloads to your Downloads folder.
5. Double-click that file and follow the installer. Click **Continue**, **Agree**, **Install**. Enter your Mac password when asked.
6. When the installer finishes, open a **new** terminal window and type:
   ```
   python3.12 --version
   ```
   You should see `Python 3.12.x`. If you do, Python is installed correctly.

### Install Python 3.12 on Windows

1. Open a web browser and go to: **https://www.python.org/downloads/**
2. Scroll down to find **Python 3.12** (look for "Looking for a specific release?" and click any `3.12.x` release).
3. Scroll down to the "Files" section and click **Windows installer (64-bit)**.
4. A file called something like `python-3.12.x-amd64.exe` downloads to your Downloads folder.
5. Double-click that file to start the installer.

   > **CRITICAL — do not skip this:** On the very first screen of the installer, there is a checkbox at the bottom that says **"Add Python 3.12 to PATH"**. **You must check this box before clicking Install.** If you miss it, uninstall Python and start over.

6. Click **Install Now** and let the installer finish.
7. Open a **new** PowerShell window and type:
   ```
   py -3.12 --version
   ```
   You should see `Python 3.12.x`.

---

## 5. Install Git

Git is the program that downloads ("clones") the pipeline code from the internet to your computer.

### Check if Git is already installed

In your terminal type:

```
git --version
```

- If you see something like `git version 2.x.x` — **you already have it, skip to Section 6.**
- If you see an error — follow the steps below.

### Install Git on a Mac

1. In your terminal, type:
   ```
   xcode-select --install
   ```
2. A pop-up window appears asking to install Command Line Developer Tools. Click **Install** and wait — this may take several minutes.
3. When it finishes, type `git --version` again to confirm it worked.

### Install Git on Windows

1. Open a web browser and go to: **https://git-scm.com/download/win**
2. The download should start automatically. If not, click the link for the 64-bit installer.
3. Run the downloaded `.exe` file. Click **Next** through every screen — the default options are all fine.
4. Open a **new** PowerShell window and type:
   ```
   git --version
   ```
   You should see a version number.

---

## 6. Download the pipeline code

"Cloning" means copying the pipeline from GitHub (a website that stores code) to your own computer.

### Where will the files go?

The pipeline will be downloaded into a folder called `Image_Analysis_Pipeline` inside whichever folder your terminal is currently in. By default that is your home folder — a fine place to put it.

### Steps (same on Mac and Windows)

1. Open a terminal.

2. Type the following and press **Enter**:
   ```
   git clone https://github.com/Theoenogo/Image_Analysis_Pipeline.git
   ```
   You will see lines of text appear as the files download. Wait for it to finish.

3. Now navigate **into** the downloaded folder. The `cd` command means "change directory" (i.e., open a folder):
   ```
   cd Image_Analysis_Pipeline
   ```

4. Confirm you are in the right place by listing the files:

   - **Mac:**
     ```
     ls
     ```
   - **Windows:**
     ```
     dir
     ```

   You should see folder names like `preprocessing`, `roi_drawing`, `roi_cropping`, `background_subtraction`, `manders_mcc`, and a file called `README.md`.

> **Tip:** Think of `cd` like double-clicking a folder to open it, except in the terminal. You can type `cd ..` at any time to go back up one level (like clicking the back button).

---

## 7. Create a virtual environment

### What is a virtual environment?

A virtual environment is a private, self-contained box for this project's Python packages. It keeps the pipeline's software completely separate from anything else on your computer, so nothing interferes with anything else. Think of it like a dedicated drawer just for this project's tools.

You only create the virtual environment **once**. After that you just "activate" it at the start of each work session (covered in Section 10).

### Make sure you are in the right folder first

Your terminal must be inside the `Image_Analysis_Pipeline` folder. If you just did Section 6, you already are. If you opened a fresh terminal, type:

```
cd Image_Analysis_Pipeline
```

(If that gives an error, you may need the full path, like `cd Documents/Image_Analysis_Pipeline` — adjust based on where you saved it.)

### Create the environment

**Mac:**
```
python3.12 -m venv venv
```

**Windows (PowerShell):**
```
py -3.12 -m venv venv
```

The command runs silently for a few seconds and then returns to the prompt. A new folder called `venv` will appear inside `Image_Analysis_Pipeline` — that is the virtual environment.

### Activate the environment

You need to activate the environment to "step inside" it before running anything.

**Mac:**
```
source venv/bin/activate
```

**Windows (PowerShell):**
```
venv\Scripts\Activate.ps1
```

After running this, you will see `(venv)` appear at the beginning of your prompt line — for example:

```
(venv) user@computer Image_Analysis_Pipeline %
```

That `(venv)` tag confirms the environment is active. **Every command from this point on should be run with `(venv)` showing.**

### Windows only — if PowerShell blocks the activation script

If Windows shows an error about "running scripts is disabled", type this once and press Enter, then try the activate command again:

```
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

Type `Y` and press Enter when prompted.

---

## 8. Install the pipeline's dependencies

### What are dependencies?

Dependencies are extra Python packages (add-ons) that the pipeline needs to work — things like image-processing libraries. The `requirements.txt` file in the project folder is a list of all of them.

### Install them (same on Mac and Windows)

Make sure `(venv)` is showing in your prompt, then type:

```
python -m pip install -r requirements.txt
```

A lot of text will scroll by as packages download and install. This can take **5–15 minutes** depending on your internet speed. That is normal — just wait.

When it finishes, the prompt will return. You will see a message like `Successfully installed ...` near the bottom.

### Verify the install worked

```
python -m pip list
```

This prints all installed packages. You should see names like `cellpose`, `numpy`, `tifffile`, `scikit-image`, and `roifile` in the list.

> **You only need to do Sections 6, 7, and 8 once.** After the initial setup, skip straight to Section 10 each time you sit down to work.

---

## 9. Running the pipeline — step by step

The pipeline has 6 steps. Steps 1–2 are run once per experiment. Step 3 is a manual check in ImageJ. Steps 4–6 can be run together in one command using the batch runner.

Before running any step, make sure:
- Your terminal shows `(venv)` at the start of the line.
- You know the full path to your data folder (see the tip below).

> **What is a "path"?** A path is the address of a file or folder on your computer. On Mac it looks like `/Users/yourname/Desktop/MyExperiment`. On Windows it looks like `C:\Users\yourname\Desktop\MyExperiment`. The easiest way to find a path:
> - **Mac:** right-click the folder in Finder, hold **Option**, and click "Copy … as Pathname".
> - **Windows:** open File Explorer, navigate to the folder, click in the address bar at the top, and copy the text that appears.

---

### Step 1 — Preprocessing

This step takes the raw TIFF images from the microscope, stacks them, corrects alignment, and sharpens them (deconvolution).

You also need **two PSF files** (Point Spread Function — special calibration images). Ask your lab manager or whoever set up the microscope where these are stored.

**Navigate to the preprocessing folder:**
```
cd preprocessing
```

**Run the script:**

Mac:
```
python preprocess.py --input-dir /path/to/your/data --gfp-psf /path/to/gfp_psf.tif --cy-psf /path/to/cy_psf.tif
```

Windows (PowerShell):
```
python preprocess.py --input-dir "C:\path\to\your\data" --gfp-psf "C:\path\to\gfp_psf.tif" --cy-psf "C:\path\to\cy_psf.tif"
```

Replace the paths in quotes with the actual locations on your computer. Progress bars will show you what is happening. When done, go back to the main folder:

```
cd ..
```

---

### Step 2 — ROI Drawing (automated cell detection)

This step automatically finds cells in your images and draws outlines around them.

**Navigate to the roi_drawing folder:**
```
cd roi_drawing
```

**Run the script:**
```
python roi_detect.py
```

The script will ask you questions interactively — type your answers and press **Enter** each time. It will ask for your data folder path, which channel to use, and a few other settings. When in doubt, press **Enter** to accept the default value (shown in square brackets).

When done, go back to the main folder:
```
cd ..
```

---

### Step 3 — Manual ROI editing in ImageJ (the only manual step)

Open Fiji/ImageJ, run the macro `Color_Merge_Automated_PreloadROIs_Adjust.ijm` (found in `roi_cropping/matlab_reference/`), and check/fix the cell outlines. Save the edited outlines to the `roi/` folder inside each sample. See your lab's ImageJ documentation for detailed instructions on this step.

---

### Steps 4–6 — Batch runner (crop + background subtract + colocalization)

Once you have finished the manual ImageJ step, you can run the remaining three steps together with a single command from the main `Image_Analysis_Pipeline` folder.

Make sure you are in the main folder (not inside any subfolder):
```
cd Image_Analysis_Pipeline
```
*(If you are already there, skip this.)*

**Run:**

Mac:
```
python run_experiment.py --input-dir /path/to/your/data
```

Windows:
```
python run_experiment.py --input-dir "C:\path\to\your\data"
```

This will automatically run:
- **Step 4** — ROI Cropping (cuts out individual cells)
- **Step 5** — Background Subtraction (removes background noise)
- **Step 6** — Colocalization analysis (produces the final CSV results)

The final results will appear as `.csv` files inside your data folder, which you can open in Excel.

---

## 10. Every time you come back: the daily checklist

Every time you open a new terminal to work on the pipeline, you need to:

1. **Open a terminal** (Section 3).

2. **Navigate to the pipeline folder:**
   ```
   cd Image_Analysis_Pipeline
   ```

3. **Activate the virtual environment:**

   Mac:
   ```
   source venv/bin/activate
   ```
   Windows:
   ```
   venv\Scripts\Activate.ps1
   ```

4. **Check that `(venv)` appears** at the start of the prompt. If it does, you are ready to run scripts.

That is it. You do not need to reinstall anything — just activate and go.

---

## 11. Troubleshooting common problems

### "command not found" or "'python3.12' is not recognized"

Python is not installed or was not added to PATH.
- **Mac:** Re-run the Python 3.12 installer from python.org.
- **Windows:** Re-run the installer and make sure you checked **"Add Python 3.12 to PATH"** on the first screen.

---

### "No such file or directory" when using a path

The path you typed does not exist or has a typo.
- Double-check the path by navigating to the folder in Finder (Mac) or File Explorer (Windows).
- On Windows, make sure you wrapped the path in quotes: `"C:\Users\..."`.
- Paths are case-sensitive on Mac — `Desktop` and `desktop` are different.

---

### The activation script is blocked on Windows

You see an error about execution policy. Run this once in PowerShell:
```
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```
Type `Y` and press Enter, then try activating again.

---

### pip install fails or a package fails to install

Try upgrading pip first, then reinstall:
```
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

---

### I see `(venv)` but scripts still fail with "ModuleNotFoundError"

The virtual environment may not have the packages installed. Run:
```
python -m pip install -r requirements.txt
```

---

### A script runs but immediately prints an error about a missing file or folder

Check that:
1. The `--input-dir` path you provided actually exists and is spelled correctly.
2. Your data folder has the expected subfolders (`gfp/`, `cy/`, etc.) — ask a labmate who has run the pipeline before.

---

### I accidentally closed the terminal mid-run

Nothing is broken. Open a new terminal, navigate to the project folder, activate `(venv)`, and re-run the step that was interrupted.

---

### I am not sure which folder I am currently in

Type this to print your current location:

- Mac: `pwd`
- Windows: `cd` (with nothing after it)

---

*If none of the above helps, copy the full error message and show it to your labmate or the repository owner.*
