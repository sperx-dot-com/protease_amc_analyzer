import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
import re
import os
from typing import List, Dict, Tuple, Optional

# Set plot style for scientific publication
sns.set_theme(style="whitegrid", context="talk")

class FluorescenceAssay:
    """
    A class to analyze fluorescence kinetic assays for protease activity.
    Designed for modularity to be integrated into Shiny/Streamlit apps.
    """

    def __init__(self, filepath: str, skiprows: int = 9):
        """
        Initialize the assay analyzer.

        Args:
            filepath (str): Path to the Excel file.
            skiprows (int): Number of metadata rows to skip before the header.
        """
        self.filepath = filepath
        
        # Data storage
        self.raw_data = None
        self.clean_data = None
        
        # Grouping storage
        self.column_mapping = {}
        self.std_groups = []
        self.sample_groups = []
        
        # Calibration state
        self.calibration_model = None # {slope, intercept, r_squared}
        self.std_concs = []
        self.results_df = None
        
        # Load and preprocess immediately
        self._load_data(skiprows)
        self._preprocess_columns()

    def _load_data(self, skiprows: int):
        """Internal method to load data safely."""
        try:
            # Load data, ensuring the first column is recognized as metadata if needed
            self.raw_data = pd.read_excel(self.filepath, skiprows=skiprows)
            # Remove empty columns if any
            self.raw_data = self.raw_data.dropna(axis=1, how='all')
            # print(f"✅ Data loaded successfully: {self.raw_data.shape[0]} time points found.")
        except Exception as e:
            raise ValueError(f"Error loading file: {e}")

    def _preprocess_columns(self):
        """
        Parses column headers to separate Time, Standards, and Samples.
        CLEANS DATA to remove invisible spaces and force numeric types.
        """
        df = self.raw_data.copy()
        
        # Normalize column names: strip whitespace
        df.columns = [str(c).strip() for c in df.columns]
        
        # Identify Time column
        time_cols = [c for c in df.columns if 'time' in c.lower()]
        if not time_cols:
             raise ValueError("Could not find a 'Time' column. Check skiprows or file format.")
        time_col = time_cols[0]
        
        df.rename(columns={time_col: 'Time'}, inplace=True)

        # --- CRITICAL FIXES FOR DATA CLEANING ---
        # 1. Force every single column to be numeric.
        # errors='coerce' turns " " (spaces) or text into NaN
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
            
        # 2. Drop columns that are NOW completely empty (all NaNs)
        # This handles the trailing commas or empty columns often found in machine exports
        df = df.dropna(axis=1, how='all')
        # ----------------------------------------
        
        # Map raw column names to grouping names
        column_mapping = {}
        for col in df.columns:
            if col == 'Time' or col == 'Reading':
                continue
            
            # Regex to capture name before the well ID in parenthesis
            match = re.match(r"^(.*)\s\([A-H]\d{2}\)$", col)
            if match:
                clean_name = match.group(1).strip()
            else:
                clean_name = col 
            
            column_mapping[col] = clean_name

        self.clean_data = df
        self.column_mapping = column_mapping
        
        # Separate Standards and Samples
        unique_groups = sorted(list(set(column_mapping.values())))
        self.std_groups = [g for g in unique_groups if g.lower().startswith('std')]
        self.sample_groups = [g for g in unique_groups if not g.lower().startswith('std')]

    def run_calibration(self, start_conc: float, increment: float):
        """
        Programmatic calibration method (for Apps/APIs).
        """
        self.std_concs = [start_conc + (i * increment) for i in range(len(self.std_groups))]
        self._perform_calibration()

    def run_calibration_interactive(self):
        """
        Interactive CLI method to ask user for standard concentrations.
        """
        print("\n--- Calibration Setup ---")
        print(f"There are {len(self.std_groups)} standards in your file.")
        
        try:
            start_conc = float(input(f"Enter concentration of the first standard ({self.std_groups[0]}) [in pmol]: "))
            increment = float(input("Enter concentration increment [in pmol]: "))
            self.run_calibration(start_conc, increment)
            
        except ValueError:
            print("❌ Invalid input. Please enter numbers only.")

    def _perform_calibration(self):
        """Calculates the standard curve using the mean signal of standards."""
        std_signals = []
        
        for group in self.std_groups:
            cols = [k for k, v in self.column_mapping.items() if v == group]
            # Average across replicates (axis=1) then across time (mean)
            avg_signal = self.clean_data[cols].mean(axis=1).mean()
            std_signals.append(avg_signal)
            
        # Safety: Drop NaNs from calibration data if any standard failed completely
        # This prevents NaN propagation to the model
        valid_indices = ~np.isnan(std_signals)
        clean_concs = np.array(self.std_concs)[valid_indices]
        clean_signals = np.array(std_signals)[valid_indices]

        if len(clean_concs) < 2:
             raise ValueError("Not enough valid standard data points for calibration.")

        slope, intercept, r_value, p_value, std_err = linregress(clean_concs, clean_signals)
        
        self.calibration_model = {
            "slope": slope,
            "intercept": intercept,
            "r_squared": r_value**2,
            "equation": f"RFU = {slope:.2f} * [AMC] + {intercept:.2f}"
        }

    def analyze_kinetics(self, start_time: float, end_time: float):
        """Analyzes sample activity within a specific time window."""
        if not self.calibration_model:
            raise RuntimeError("Calibration not found. Run run_calibration() first.")
    
        mask = (self.clean_data['Time'] >= start_time) & (self.clean_data['Time'] <= end_time)
        df_window = self.clean_data.loc[mask]
        
        results = []

        for group in self.sample_groups:
            cols = [k for k, v in self.column_mapping.items() if v == group]
            # The previous error happened here because 'mean_trace' had strings
            mean_trace = df_window[cols].mean(axis=1)
            times = df_window['Time']
            
            # --- ROBUST REGRESSION ---
            # Remove any NaNs from this specific trace before regression
            # (e.g., if one timepoint is bad)
            valid_idx = ~np.isnan(mean_trace)
            
            if np.sum(valid_idx) > 1:
                slope_rfu, intercept_rfu, r_sq, _, _ = linregress(times[valid_idx], mean_trace[valid_idx])
            else:
                # Fallback if trace is empty or all NaNs
                slope_rfu, intercept_rfu, r_sq = 0.0, 0.0, 0.0
            
            rate_pmol_sec = slope_rfu / self.calibration_model['slope']
            
            # Use last valid point for endpoint
            if np.sum(valid_idx) > 0:
                end_rfu = mean_trace[valid_idx].iloc[-1]
                amount_pmol = (end_rfu - self.calibration_model['intercept']) / self.calibration_model['slope']
            else:
                amount_pmol = 0.0
            
            results.append({
                "Sample": group,
                "Rate_RFU_s": slope_rfu,
                "Rate_pmol_s": rate_pmol_sec,
                "Endpoint_pmol": amount_pmol,
                "Linearity_R2": r_sq
            })
            
        self.results_df = pd.DataFrame(results)
        # Final cleanup: Replace any remaining NaNs in results with 0 or None so JSON serialization works
        self.results_df = self.results_df.fillna(0)

    def export_excel(self, filename="analysis_results.xlsx"):
        """Export results and calibration data to an Excel file."""
        if self.results_df is None:
            return
        with pd.ExcelWriter(filename) as writer:
            self.results_df.to_excel(writer, sheet_name='Results', index=False)
            
            calib_data = {
                "Parameter": ["Slope", "Intercept", "R_Squared", "Equation"],
                "Value": [
                    self.calibration_model['slope'],
                    self.calibration_model['intercept'],
                    self.calibration_model['r_squared'],
                    self.calibration_model['equation']
                ]
            }
            pd.DataFrame(calib_data).to_excel(writer, sheet_name='Calibration', index=False)

    def save_individual_plots(self, output_dir="plots"):
        """Generates and saves individual PNG files for reports."""
        if self.results_df is None:
            return

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        # 1. Standard Curve
        fig1 = plt.figure(figsize=(8, 6))
        x_pred = np.linspace(min(self.std_concs), max(self.std_concs), 100)
        y_pred = self.calibration_model['slope'] * x_pred + self.calibration_model['intercept']
        std_signals = [self.clean_data[[k for k,v in self.column_mapping.items() if v == g]].mean().mean() for g in self.std_groups]
        
        plt.scatter(self.std_concs, std_signals, color='blue', s=100, label='Standards')
        plt.plot(x_pred, y_pred, 'r--', label=f"Fit ($R^2={self.calibration_model['r_squared']:.3f}$)")
        plt.title("Standard Curve (AMC)")
        plt.xlabel("Concentration [pmol]")
        plt.ylabel("RFU")
        plt.legend()
        plt.tight_layout()
        fig1.savefig(f"{output_dir}/standard_curve.png", dpi=300)
        plt.close(fig1)
        
        # 2. Kinetic Traces
        fig2 = plt.figure(figsize=(10, 6))
        for group in self.sample_groups:
            cols = [k for k, v in self.column_mapping.items() if v == group]
            mean_trace = self.clean_data[cols].mean(axis=1)
            plt.plot(self.clean_data['Time'], mean_trace, label=group)
        plt.title("Kinetic Traces (Averaged)")
        plt.xlabel("Time [s]")
        plt.ylabel("RFU")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        fig2.savefig(f"{output_dir}/kinetic_traces.png", dpi=300)
        plt.close(fig2)
        
        # 3. Activity Bar Plot
        fig3 = plt.figure(figsize=(8, 6))
        sns.barplot(data=self.results_df, x="Sample", y="Rate_pmol_s", palette="viridis", hue="Sample", legend=False)
        plt.title("Protease Activity")
        plt.ylabel("Rate [pmol AMC / sec]")
        plt.xticks(rotation=45)
        plt.tight_layout()
        fig3.savefig(f"{output_dir}/activity_barplot.png", dpi=300)
        plt.close(fig3)
        

    def plot_overview(self):
        """
        Generates a composite figure of the analysis.
        Returns the matplotlib Figure object (for use in Apps).
        """
        if self.results_df is None:
            return None

        # Increase figure size slightly for better dashboard look
        fig = plt.figure(figsize=(16, 10))
        gs = fig.add_gridspec(2, 2)

        # 1. Standard Curve
        ax1 = fig.add_subplot(gs[0, 0])
        x_pred = np.linspace(min(self.std_concs), max(self.std_concs), 100)
        y_pred = self.calibration_model['slope'] * x_pred + self.calibration_model['intercept']
        
        std_signals = [self.clean_data[[k for k,v in self.column_mapping.items() if v == g]].mean().mean() for g in self.std_groups]
        
        ax1.scatter(self.std_concs, std_signals, color='blue', s=100, label='Standards')
        ax1.plot(x_pred, y_pred, 'r--', label=f"Fit ($R^2={self.calibration_model['r_squared']:.4f}$)")
        ax1.set_title("Standard Curve (AMC)")
        ax1.set_xlabel("Concentration [pmol]")
        ax1.set_ylabel("RFU")
        ax1.legend()

        # 2. Kinetic Traces
        ax2 = fig.add_subplot(gs[0, 1])
        for group in self.sample_groups:
            cols = [k for k, v in self.column_mapping.items() if v == group]
            mean_trace = self.clean_data[cols].mean(axis=1)
            ax2.plot(self.clean_data['Time'], mean_trace, label=group)
        
        ax2.set_title("Kinetic Traces (Averaged Replicates)")
        ax2.set_xlabel("Time [s]")
        ax2.set_ylabel("RFU")
        ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

        # 3. Calculated Activity Bar Chart
        ax3 = fig.add_subplot(gs[1, 0])
        sns.barplot(
            data=self.results_df, 
            x="Sample", 
            y="Rate_pmol_s", 
            ax=ax3, 
            palette="viridis", 
            hue="Sample", 
            legend=False
        )
        ax3.set_title("Protease Activity")
        ax3.set_ylabel("Rate [pmol AMC / sec]")
        ax3.tick_params(axis='x', rotation=45)

        # 4. Results Table
        ax4 = fig.add_subplot(gs[1, 1])
        ax4.axis('tight')
        ax4.axis('off')
        
        table_df = self.results_df[['Sample', 'Rate_pmol_s', 'Endpoint_pmol', 'Linearity_R2']].copy()
        table_df['Rate_pmol_s'] = table_df['Rate_pmol_s'].round(3)
        table_df['Endpoint_pmol'] = table_df['Endpoint_pmol'].round(2)
        table_df['Linearity_R2'] = table_df['Linearity_R2'].round(4)
        
        table_data = table_df.values.tolist()
        col_labels = ["Sample", "Rate\n(pmol/s)", "Endpoint\n(pmol)", "$R^2$"]
        
        table = ax4.table(cellText=table_data, colLabels=col_labels, loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5)
        ax4.set_title("Results Summary", pad=20)
        
        plt.tight_layout()
        return fig

# --- Main Execution Block (CLI Mode) ---
if __name__ == "__main__":
    assay = None
    
    # 1. Initialize with Drag & Drop support
    print("\n--- Protease Analysis Tool ---")
    while True:
        filepath = input("Enter path to Excel file (or drag & drop file) [or 'q' to quit]: ").strip()
        
        # Handle quit
        if filepath.lower() in ['q', 'quit', 'exit']:
            print("Exiting...")
            exit()
            
        # Remove quotes if present (common in drag & drop on Windows/Mac)
        filepath = filepath.strip('"').strip("'")
        
        if not filepath:
            continue

        try:
            print(f"\nLoading data from: {filepath}...")
            assay = FluorescenceAssay(filepath)
            print("✅ File loaded successfully.")
            break
            
        except FileNotFoundError:
            print(f"❌ Error: File not found at '{filepath}'")
            print("   Please check the path and try again.")
        except ValueError as e:
            print(f"❌ Data Error: {e}")
            print("   Ensure this is the correct Excel file format.")
        except Exception as e:
            print(f"❌ Unexpected error: {e}")
            
    # 2. Run Calibration
    if assay:
        try:
            assay.run_calibration_interactive()
            
            # 3. Run Analysis
            print("\n--- Analysis Parameters ---")
            t_start = float(input("Enter start time for linear range [s] (e.g., 0): "))
            t_end = float(input("Enter end time for linear range [s] (e.g., 200): "))
            
            assay.analyze_kinetics(start_time=t_start, end_time=t_end)
            
            # 4. Export Options
            print("\n--- Export Options ---")
            
            # Robust Excel Export Loop
            while True:
                try:
                    choice = input("Export results to Excel? (y/n): ").lower().strip()
                    if choice != 'y':
                        break
                    
                    out_file = input("Enter filename (default: results.xlsx): ").strip()
                    if not out_file: 
                        out_file = "results.xlsx"
                    
                    if not out_file.lower().endswith('.xlsx'):
                        out_file += '.xlsx'
                        
                    assay.export_excel(out_file)
                    print(f"✅ Successfully exported to {out_file}")
                    break
                    
                except Exception as e:
                    print(f"❌ Error exporting file: {e}")
                    print("   Please try again.")

            # Robust Plot Export Loop
            while True:
                try:
                    choice = input("Save individual plots as PNGs? (y/n): ").lower().strip()
                    if choice != 'y':
                        break
                    
                    out_folder = input("Enter folder name (default: plots): ").strip()
                    if not out_folder: 
                        out_folder = "plots"
                    
                    assay.save_individual_plots(out_folder)
                    print(f"✅ Plots saved to folder '{out_folder}/'")
                    break
                    
                except Exception as e:
                    print(f"❌ Error saving plots: {e}")
                    print("   Please try again.")
            
            # 5. Visualize
            print("\nGenerating overview plot...")
            fig = assay.plot_overview()
            if fig:
                plt.show()

        except KeyboardInterrupt:
            print("\nAnalysis cancelled by user.")
        except Exception as e:
            print(f"\n❌ An error occurred during analysis: {e}")