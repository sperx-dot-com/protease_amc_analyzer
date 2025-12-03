from shiny import App, ui, render, reactive
import pandas as pd
import numpy as np
from protease_amc_analyzer import FluorescenceAssay

import plotly.graph_objects as go
from plotly.subplots import make_subplots
from shinywidgets import output_widget, render_widget
# --- UI Definition ---
app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.h4("Experiment Settings"),
        
        ui.input_file("file_upload", "Upload Excel File (.xlsx)", accept=[".xlsx", ".xls"]),
        
        ui.hr(),
        ui.h5("Calibration"),
        ui.input_numeric("cal_start", "Std Start (pmol)", value=0),
        ui.input_numeric("cal_inc", "Std Increment (pmol)", value=20),
        
        ui.hr(),
        ui.h5("Kinetics"),
        ui.input_numeric("time_start", "Start Time (s)", value=0),
        ui.input_numeric("time_end", "End Time (s)", value=200),
        
        ui.hr(),
        ui.h5("Meta"),
        ui.input_numeric("skip_rows", "Header Rows to Skip", value=9),
    ),
    
    ui.layout_columns(
        ui.card(
            ui.card_header("Interactive Analysis Overview"),
            # output_widget replaces output_plot for interactive graphs
            output_widget("plot_overview", height="900px")
        ),
        col_widths=12
    ),
    
    ui.layout_columns(
        ui.card(
            ui.card_header("Calculated Results"),
            ui.output_data_frame("results_table")
        ),
        ui.card(
            ui.card_header("Calibration Stats"),
            ui.output_text_verbatim("calibration_stats")
        )
    ),
    fillable=False
)

# --- Server Logic ---
def server(input, output, session):
    
    @reactive.calc
    def get_assay():
        file_info = input.file_upload()
        if not file_info: return None
        try:
            return FluorescenceAssay(file_info[0]["datapath"], skiprows=input.skip_rows())
        except Exception: return None

    @reactive.calc
    def run_analysis():
        assay = get_assay()
        if not assay: return None
        
        assay.run_calibration(start_conc=input.cal_start(), increment=input.cal_inc())
        assay.analyze_kinetics(start_time=input.time_start(), end_time=input.time_end())
        return assay

    # --- Interactive Plotting Logic ---
    @render_widget
    def plot_overview():
        assay = run_analysis()
        if not assay or assay.results_df is None:
            return go.Figure() # Return empty figure if no data
        
        # Create a 2x2 subplot layout
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=("Standard Curve", "Kinetic Traces", "Protease Activity", "Results Table"),
            specs=[[{"type": "xy"}, {"type": "xy"}],
                   [{"type": "xy"}, {"type": "table"}]], # 4th is a table type
            vertical_spacing=0.15
        )

        # 1. Standard Curve (Scatter + Line)
        # Get actual signals
        std_signals = [assay.clean_data[[k for k,v in assay.column_mapping.items() if v == g]].mean().mean() for g in assay.std_groups]
        
        # Points
        fig.add_trace(go.Scatter(
            x=assay.std_concs, y=std_signals,
            mode='markers', name='Standards', marker=dict(color='blue', size=10)
        ), row=1, col=1)
        
        # Fit Line
        x_pred = np.linspace(min(assay.std_concs), max(assay.std_concs), 100)
        y_pred = assay.calibration_model['slope'] * x_pred + assay.calibration_model['intercept']
        fig.add_trace(go.Scatter(
            x=x_pred, y=y_pred,
            mode='lines', name=f'Fit (R²={assay.calibration_model["r_squared"]:.4f})',
            line=dict(color='red', dash='dash')
        ), row=1, col=1)

        # 2. Kinetic Traces (Lines)
        # We loop through samples and add a trace for each
        for group in assay.sample_groups:
            cols = [k for k, v in assay.column_mapping.items() if v == group]
            mean_trace = assay.clean_data[cols].mean(axis=1)
            fig.add_trace(go.Scatter(
                x=assay.clean_data['Time'], y=mean_trace,
                mode='lines', name=group
            ), row=1, col=2)

        # 3. Activity Bar Chart
        fig.add_trace(go.Bar(
            x=assay.results_df['Sample'],
            y=assay.results_df['Rate_pmol_s'],
            name='Activity', marker=dict(color=assay.results_df['Rate_pmol_s'], colorscale='Viridis')
        ), row=2, col=1)

        # 4. Results Table (embedded in plot)
        df_display = assay.results_df.copy()
        fig.add_trace(go.Table(
            header=dict(values=["Sample", "Rate (pmol/s)", "End (pmol)", "R²"],
                        fill_color='paleturquoise', align='left'),
            cells=dict(values=[df_display.Sample, 
                               df_display.Rate_pmol_s.round(3), 
                               df_display.Endpoint_pmol.round(1), 
                               df_display.Linearity_R2.round(4)],
                       fill_color='lavender', align='left')
        ), row=2, col=2)

        # Update Layout for a clean look
        fig.update_layout(
            height=900, 
            showlegend=True,
            title_text="Protease Assay Analysis Dashboard",
            template="plotly_white"
        )
        
        # Update axes labels
        fig.update_xaxes(title_text="Concentration [pmol]", row=1, col=1)
        fig.update_yaxes(title_text="RFU", row=1, col=1)
        fig.update_xaxes(title_text="Time [s]", row=1, col=2)
        fig.update_yaxes(title_text="RFU", row=1, col=2)
        fig.update_yaxes(title_text="Rate [pmol/s]", row=2, col=1)

        return fig

    @render.data_frame
    def results_table():
        assay = run_analysis()
        if not assay or assay.results_df is None: return pd.DataFrame()
        
        # Return nice table
        df_display = assay.results_df.copy()
        df_display['Rate_pmol_s'] = df_display['Rate_pmol_s'].round(3)
        df_display['Endpoint_pmol'] = df_display['Endpoint_pmol'].round(1)
        df_display['Linearity_R2'] = df_display['Linearity_R2'].round(4)
        cols = ["Sample", "Rate_pmol_s", "Endpoint_pmol", "Linearity_R2"]
        return render.DataGrid(df_display[cols], width="100%")

    @render.text
    def calibration_stats():
        assay = run_analysis()
        if not assay or not assay.calibration_model: return "Upload file to see stats."
        m = assay.calibration_model
        return f"Equation: {m['equation']}\nR²: {m['r_squared']:.4f}"

app = App(app_ui, server)

if __name__ == "__main__":
    app.run()