# **Protease Kinetics Analyzer**

A robust, object-oriented Python tool for analyzing fluorescence kinetic assays, specifically designed for protease activity measurements (e.g., AMC release) using data exported from **ThermoScientific Varioskan Lux**.  
This project provides both a **Command Line Interface (CLI)** for quick automated scripts and a modern **Shiny for Python Web App** for interactive data exploration.

## **🧪 Scientific Context**

This tool is designed for **enzyme kinetics** workflows:

1. **Input:** Time-resolved fluorescence data (.xlsx) exported from SkanIt software.  
2. **Calibration:** automatically generates a standard curve from Std wells (Linear Regression: $RFU \= m \\cdot \[Conc\] \+ c$).  
3. **Kinetics:** Calculates the initial reaction velocity ($v\_0$) and endpoint product formation for unknown samples.  
4. **Output:** Calculated rates (pmol/s), $R^2$ linearity checks, and publication-quality plots.

## **📦 Installation**

Ensure you have Python 3.8+ installed. Install the required dependencies:  
pip install pandas numpy matplotlib seaborn scipy shiny shinywidgets plotly openpyxl

## **🚀 How to Run**

### **1\. Interactive Web Dashboard (Recommended)**

Launch the Shiny app to drag-and-drop files, adjust time windows interactively, and explore data with zoomable Plotly graphs.  
python app.py

*Then open the local URL displayed in your terminal (usually http://127.0.0.1:8000).*

### **2\. Command Line Interface (CLI)**

Run the analysis directly in the terminal. The script will guide you through loading the file and setting parameters.  
python protease\_analysis.py

### **3\. Use as a Python Library**

You can import the core logic into Jupyter Notebooks or other scripts:  
from protease\_analysis import FluorescenceAssay

\# 1\. Initialize  
assay \= FluorescenceAssay("my\_data.xlsx", skiprows=9)

\# 2\. Calibrate (e.g., 0 to 100 pmol)  
assay.run\_calibration(start\_conc=0, increment=20)

\# 3\. Analyze Kinetics (e.g., 0 to 200 seconds)  
assay.analyze\_kinetics(start\_time=0, end\_time=200)

\# 4\. Export  
assay.export\_excel("results.xlsx")

## **📂 Project Structure**

The project follows strict **Separation of Concerns**:

* **protease\_analysis.py (The Logic):**  
  * Contains the FluorescenceAssay class.  
  * Handles data parsing, cleaning (removing NaNs/spaces), math, and statistics.  
  * **Zero UI dependencies:** Can run on a server or headless environment.  
  * Includes a robust cleaning pipeline to handle ThermoScientific export quirks (trailing commas, empty columns).  
* **app.py (The Interface):**  
  * A **Shiny for Python** application.  
  * Handles user interaction (Sliders, File Uploads).  
  * Uses **Plotly** for interactive visualization.  
  * Connects to protease\_analysis.py to perform calculations.

## **📊 Input File Format**

The tool expects standard **Excel exports (.xlsx)** from SkanIt software with the following structure:

* **Header Rows:** \~9 rows of metadata (instrument settings) to be skipped.  
* **Time Column:** A column containing time (seconds).  
* **Data Columns:**  
  * **Standards:** Must be named starting with Std (e.g., Std0001, Std0002). Replicates (e.g., Std0001 (A01)) are automatically grouped and averaged.  
  * **Samples:** Any other column name. Replicates with the same name are automatically grouped.

## **🛡️ robust Features**

* **Nan-Safe Regression:** Automatically ignores individual bad data points (NaN or empty cells) within a kinetic trace without crashing the analysis.  
* **Auto-Cleaning:** Detects and removes empty columns or text-formatted numbers often found in raw machine exports.

*Created for biotechnological data science workflows.*
