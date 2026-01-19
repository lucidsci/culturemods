"""
NiceGUI application for configuring, running, and visualizing
1D oxygen diffusion-reaction simulations.
"""

from nicegui import ui, run
import asyncio
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from rxd_fipy_1d import SimulationConfig, SimulationResult, run_simulation
from conversions import rate_pmols_per_L_per_minute_to_umolar_per_s, flux_fmols_per_mm2_per_s_to_umolar_per_s

class SimulationGUI:
    """GUI for oxygen diffusion simulation."""

    def __init__(self):
        self.result: SimulationResult | None = None
        self.df = None

        # Default values
        self.config_values = {
            'D': 3e-3,
            'C_air': 200.0,
            'L': 3.1,
            'nz': 100,
            'rate': 5, #pmols/L/min
            'k1': 0.0,
            'flux_bottom': 0.0,  # fmol/mm²/s
            'top_constraint': 'open',
            'dt': 1.0,
            'steps': 1800,
            'C_initial_fraction': 1.0,
            'halt_on_C_zero': True,
        }

        # Visualization parameters
        self.probe_height_mm = 1.0
        self.profile_step = 0

        # UI element references
        self.timeseries_plot = None
        self.profile_plot = None
        self.status_label = None
        self.probe_slider = None
        self.step_slider = None
        self.progress_bar = None
        self.run_button = None

        # Progress tracking
        self.current_step = 0
        self.total_steps = 0
        self.is_running = False

    def build_ui(self):
        """Build the main UI layout."""
        ui.dark_mode().enable()

        with ui.header().classes('items-center justify-between'):
            ui.label('Culture Well O2 Diffusion-Reaction Simulator').classes('text-2xl font-bold')

        with ui.row().classes('w-full gap-4 p-4'):
            # Left panel - Configuration
            with ui.card().classes('w-80'):
                ui.label('Simulation Parameters').classes('text-lg font-bold mb-2')
                self._build_config_panel()

            # Right panel - Visualization
            with ui.column().classes('flex-grow gap-4'):
                self._build_visualization_panel()

    def _build_config_panel(self):
        """Build the configuration input panel."""
        with ui.column().classes('gap-2 w-full'):
            # Physical parameters
            ui.label('Physical Parameters').classes('font-semibold text-blue-400')

            ui.number(
                'Diffusion coeff D (mm²/s)',
                value=self.config_values['D'],
                format='%.2e',
                step=1e-4,
                on_change=lambda e: self._update_config('D', e.value)
            ).classes('w-full')

            ui.number(
                'C_air (µM)',
                value=self.config_values['C_air'],
                format='%.1f',
                on_change=lambda e: self._update_config('C_air', e.value)
            ).classes('w-full')

            ui.number(
                'Well depth L (mm)',
                value=self.config_values['L'],
                format='%.2f',
                step=0.1,
                on_change=lambda e: self._update_config('L', e.value)
            ).classes('w-full')

            ui.number(
                'Mesh points (nz)',
                value=self.config_values['nz'],
                format='%.0f',
                min=10,
                max=500,
                step=10,
                on_change=lambda e: self._update_config('nz', int(e.value))
            ).classes('w-full')

            ui.separator()

            # Reaction parameters
            ui.label('Reaction Parameters').classes('font-semibold text-blue-400')

            ui.number(
                'Zero-order reaction rate (pmols/L/minute))',
                value=self.config_values['rate'],
                format='%.1f',
                step=.1,
                on_change=lambda e: self._update_config('rate', e.value)
            ).classes('w-full')

            ui.number(
                'First-order rate k1',
                value=self.config_values['k1'],
                format='%.6f',
                step=0.001,
                min=0,
                on_change=lambda e: self._update_config('k1', e.value)
            ).classes('w-full')

            ui.number(
                'Bottom flux (fmol/mm²/s)',
                value=self.config_values['flux_bottom'],
                format='%.2f',
                step=0.1,
                min=0,
                on_change=lambda e: self._update_config('flux_bottom', e.value)
            ).classes('w-full')

            ui.select(
                ['open', 'sealed'],
                value=self.config_values['top_constraint'],
                label='Top boundary',
                on_change=lambda e: self._update_config('top_constraint', e.value)
            ).classes('w-full')

            ui.separator()

            # Time parameters
            ui.label('Time Parameters').classes('font-semibold text-blue-400')

            ui.number(
                'Time step dt (s)',
                value=self.config_values['dt'],
                format='%.2f',
                step=0.1,
                min=0.01,
                on_change=lambda e: self._update_config('dt', e.value)
            ).classes('w-full')

            ui.number(
                'Number of steps',
                value=self.config_values['steps'],
                format='%.0f',
                min=1,
                max=10000,
                step=10,
                on_change=lambda e: self._update_config('steps', int(e.value))
            ).classes('w-full')

            ui.number(
                'Initial C fraction',
                value=self.config_values['C_initial_fraction'],
                format='%.2f',
                min=0,
                max=1,
                step=0.1,
                on_change=lambda e: self._update_config('C_initial_fraction', e.value)
            ).classes('w-full')

            ui.checkbox(
                'Halt when C=0',
                value=self.config_values['halt_on_C_zero'],
                on_change=lambda e: self._update_config('halt_on_C_zero', e.value)
            )

            ui.separator()

            # Run button and status
            with ui.row().classes('w-full gap-2'):
                self.run_button = ui.button('Run Simulation', on_click=self._run_simulation).classes('flex-grow')

            self.progress_bar = ui.linear_progress(value=0, show_value=False).classes('w-full')
            self.progress_bar.visible = False

            self.status_label = ui.label('Ready').classes('text-sm text-gray-400')

            # Derived parameters display
            with ui.expansion('Derived Parameters', icon='info').classes('w-full'):
                self.derived_params_container = ui.column().classes('gap-1')
                self._update_derived_display()

    def _build_visualization_panel(self):
        """Build the visualization panel with plots."""
        # Timeseries plot card
        with ui.card().classes('w-full'):
            with ui.row().classes('items-center gap-4 mb-2'):
                ui.label('O₂ Concentration vs Time').classes('text-lg font-bold')
                ui.label('Probe height:').classes('text-sm')
                self.probe_slider = ui.slider(
                    min=0, max=3.0, step=0.1, value=self.probe_height_mm,
                    on_change=self._on_probe_height_change
                ).classes('w-48')
                self.probe_label = ui.label(f'{self.probe_height_mm:.1f} mm').classes('text-sm w-16')

            self.timeseries_plot = ui.plotly({}).classes('w-full h-[400px]')

        # Profile heatmap card
        with ui.card().classes('w-full'):
            with ui.row().classes('items-center gap-4 mb-2'):
                ui.label('Vertical Concentration Profile').classes('text-lg font-bold')
                ui.label('Time step:').classes('text-sm')
                self.step_slider = ui.slider(
                    min=0, max=self.config_values['steps'], step=1, value=self.profile_step,
                    on_change=self._on_step_change
                ).classes('w-48')
                self.step_label = ui.label('Step 0 (0.0 s)').classes('text-sm w-32')

            with ui.row().classes('w-full gap-4'):
                # 1D heatmap
                self.profile_plot = ui.plotly({}).classes('flex-grow h-[400px]')

        # Initialize empty plots
        self._update_timeseries_plot()
        self._update_profile_plot()

    def _update_config(self, key: str, value):
        """Update configuration value."""
        self.config_values[key] = value

        # Update probe slider max based on well depth
        if key == 'L' and self.probe_slider:
            self.probe_slider.props['max']= value
            if self.probe_height_mm > value:
                self.probe_height_mm = value
                self.probe_slider.value = value

        print(f'config update {key} = {value}')
        # Update step slider max based on steps
        if key == 'steps' and self.step_slider:
            self.step_slider.props['max'] = value

        self._update_derived_display()

    def _update_derived_display(self):
        """Update the derived parameters display."""
        if not hasattr(self, 'derived_params_container'):
            return

        self.derived_params_container.clear()
        with self.derived_params_container:
            try:
                config = self._create_config()
                ui.label(f'dz: {config.dz:.4f} mm').classes('text-xs')
                ui.label(f'T (char. time): {config.T:.1f} s').classes('text-xs')
                ui.label(f'Total time: {config.total_time:.1f} s ({config.total_time/3600:.2f} hr)').classes('text-xs')
                ui.label(f'Damköhler #: {config.damkohler:.3f}').classes('text-xs')
                if config.L_char != float('inf'):
                    ui.label(f'L_char: {config.L_char:.3f} mm').classes('text-xs')
            except Exception:
                ui.label('Invalid parameters').classes('text-xs text-red-400')

    def _create_config(self) -> SimulationConfig:
        """Create SimulationConfig from current values."""
        C_air = self.config_values['C_air']
        L = self.config_values['L']

        # Convert volumetric reaction rate
        rate_umolar_per_s = rate_pmols_per_L_per_minute_to_umolar_per_s(self.config_values['rate'])
        k = rate_umolar_per_s / C_air  # normalize out concentration units

        dt = self.config_values['dt']
        # Convert bottom flux (fmol/mm2/s to umols/mm2/s)
        #FIXME - don't think this is correct
        flux_bottom = self.config_values['flux_bottom'] / C_air * dt


        return SimulationConfig(
            D=self.config_values['D'],
            C_air=self.config_values['C_air'],
            L=L,
            nz=self.config_values['nz'],
            k=k,
            k1=self.config_values['k1'],
            flux_bottom=flux_bottom,
            top_constraint=self.config_values['top_constraint'],
            dt=self.config_values['dt'],
            steps=self.config_values['steps'],
            C_initial_fraction=self.config_values['C_initial_fraction'],
            halt_on_C_zero=self.config_values['halt_on_C_zero'],
        )

    def _step_callback(self, step: int):
        """Callback called on each simulation step to update progress."""
        self.current_step = step

    async def _run_simulation(self):
        """Run the simulation with current configuration."""
        if self.is_running:
            return

        self.is_running = True
        self.run_button.disable()
        self.status_label.text = 'Running simulation...'
        self.status_label.classes('text-yellow-400', remove='text-gray-400 text-green-400 text-red-400')

        # Show and reset progress bar
        self.progress_bar.visible = True
        self.progress_bar.value = 0
        self.current_step = 0

        try:
            config = self._create_config()
            self.total_steps = config.steps

            # Start progress update timer
            progress_timer = ui.timer(0.1, self._update_progress)

            # Run simulation in background thread to not block UI
            self.result = await run.io_bound(
                run_simulation,
                config,
                record_every=1,
                verbose=False,
                step_callback=self._step_callback
            )

            # Stop progress timer and set to 100%
            progress_timer.cancel()
            self.progress_bar.value = 1.0

            self.df = self.result.to_dataframe()

            # Update slider ranges
            if len(self.df) > 0:
                max_step = self.df['step'].max()
                self.step_slider.max = max_step
                if self.profile_step > max_step:
                    self.profile_step = max_step
                    self.step_slider.value = max_step
                self.probe_slider.max = config.L

            # Update plots
            self._update_timeseries_plot()
            self._update_profile_plot()

            actual_steps = len(self.df['step'].unique()) if len(self.df) > 0 else 0
            self.status_label.text = f'Completed ({actual_steps} steps recorded)'
            self.status_label.classes('text-green-400', remove='text-gray-400 text-yellow-400 text-red-400')

        except Exception as e:
            self.status_label.text = f'Error: {str(e)}'
            self.status_label.classes('text-red-400', remove='text-gray-400 text-yellow-400 text-green-400')

        finally:
            self.is_running = False
            self.run_button.enable()
            # Hide progress bar after a short delay
            await asyncio.sleep(1.0)
            self.progress_bar.visible = False

    def _update_progress(self):
        """Update the progress bar value."""
        if self.total_steps > 0:
            self.progress_bar.value = self.current_step / self.total_steps

    def _on_probe_height_change(self, e):
        """Handle probe height slider change."""
        self.probe_height_mm = e.value
        self.probe_label.text = f'{self.probe_height_mm:.1f} mm'
        self._update_timeseries_plot()

    def _on_step_change(self, e):
        """Handle time step slider change."""
        self.profile_step = int(e.value)
        if self.result:
            t_s = self.profile_step * self.result.config.dt
            self.step_label.text = f'Step {self.profile_step} ({t_s:.1f} s)'
        self._update_profile_plot()

    def _update_timeseries_plot(self):
        """Update the timeseries plot."""
        if self.df is None or len(self.df) == 0:
            # Empty plot
            fig = go.Figure()
            t_max = self.config_values['steps'] * self.config_values['dt']
            fig.update_layout(
                template='plotly_dark',
                margin=dict(l=60, r=20, t=30, b=50),
                xaxis_title='Time (minutes)',
                yaxis_title='O₂ Concentration (µM)',
                xaxis=dict(range=[0, t_max]),
                yaxis=dict(range=[0, 220]),
            )
            fig.add_annotation(
                text="Run simulation to see results",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=14, color="gray")
            )
            self.timeseries_plot.update_figure(fig)
            return

        # Find the z_idx closest to probe height
        config = self.result.config
        probe_idx = int(self.probe_height_mm / config.dz)
        probe_idx = max(0, min(probe_idx, config.nz - 1))

        df_probe = self.df[self.df['z_idx'] == probe_idx].sort_values('t_hrs')

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=df_probe['t_mins'],
            y=df_probe['C'],
            mode='lines',
            name=f'z = {self.probe_height_mm:.1f} mm',
            line=dict(color='#00d4aa', width=2)
        ))

        c_max = self.df.C.max()
        fig.update_layout(
            template='plotly_dark',
            margin=dict(l=60, r=20, t=30, b=50),
            xaxis_title='Time (minutes)',
            yaxis_title='O₂ Concentration (µM)',
            yaxis=dict(range=[0, c_max * 1.1]),
            showlegend=True,
            legend=dict(x=0.02, y=0.98),
        )

        # Add reference line for C_air
        fig.add_hline(y=config.C_air, line_dash="dash", line_color="gray",
                      annotation_text=f"C_air = {config.C_air} µM")

        self.timeseries_plot.update_figure(fig)

    def _update_profile_plot(self):
        """Update the vertical profile heatmap plot."""
        if self.df is None or len(self.df) == 0:
            # Empty plot
            fig = go.Figure()
            fig.update_layout(
                template='plotly_dark',
                margin=dict(l=60, r=20, t=30, b=50),
            )
            fig.add_annotation(
                text="Run simulation to see results",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=14, color="gray")
            )
            self.profile_plot.update_figure(fig)
            return

        config = self.result.config
        profile_df = self.result.get_profile_at_step(self.profile_step)

        if len(profile_df) == 0:
            return

        heights = profile_df['height_mm'].values
        concentrations = profile_df['C'].values

        # Create figure with two subplots: heatmap bar and line plot
        fig = make_subplots(
            rows=1, cols=2,
            column_widths=[0.15, 0.85],
            horizontal_spacing=0.02,
            shared_yaxes=True
        )

        # 1D Heatmap (vertical bar)
        # Reshape concentration for heatmap (needs 2D array)
        z_data = concentrations.reshape(-1, 1)

        fig.add_trace(go.Heatmap(
            z=z_data,
            x=[0],
            y=heights,
            colorscale='RdYlBu_r',
            zmin=0,
            zmax=config.C_air,
            showscale=False,
            hovertemplate='Height: %{y:.2f} mm<br>C: %{z:.1f} µM<extra></extra>'
        ), row=1, col=1)

        # Line plot of profile
        fig.add_trace(go.Scatter(
            x=concentrations,
            y=heights,
            mode='lines',
            line=dict(color='#00d4aa', width=2),
            hovertemplate='C: %{x:.1f} µM<br>Height: %{y:.2f} mm<extra></extra>'
        ), row=1, col=2)

        # Add marker for probe height
        probe_idx = int(self.probe_height_mm / config.dz)
        probe_idx = max(0, min(probe_idx, len(concentrations) - 1))
        if probe_idx < len(concentrations):
            fig.add_trace(go.Scatter(
                x=[concentrations[probe_idx]],
                y=[heights[probe_idx]],
                mode='markers',
                marker=dict(color='yellow', size=10, symbol='diamond'),
                name='Probe',
                hovertemplate='Probe<br>C: %{x:.1f} µM<br>Height: %{y:.2f} mm<extra></extra>'
            ), row=1, col=2)

        fig.update_layout(
            template='plotly_dark',
            margin=dict(l=60, r=20, t=30, b=50),
            showlegend=False,
        )

        c_max = np.max(concentrations)
        # Update axes
        fig.update_xaxes(showticklabels=False, row=1, col=1)
        fig.update_xaxes(title_text='O₂ Concentration (µM)', range=[0, c_max * 1.1], row=1, col=2)
        fig.update_yaxes(title_text='Height (mm)', row=1, col=1)

        # Add annotations for top/bottom
        fig.add_annotation(
            x=0, y=config.L, xref='x', yref='y',
            text="Air", showarrow=False,
            font=dict(size=10), xanchor='center', yanchor='bottom'
        )
        fig.add_annotation(
            x=0, y=0, xref='x', yref='y',
            text="Bottom", showarrow=False,
            font=dict(size=10), xanchor='center', yanchor='top'
        )

        self.profile_plot.update_figure(fig)


def main():
    """Main entry point."""
    app = SimulationGUI()
    app.build_ui()
    ui.run(title='O₂ Diffusion Simulator', port=8081, reload=False)


if __name__ in {'__main__', '__mp_main__'}:
    main()
