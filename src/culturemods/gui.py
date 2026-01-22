"""
NiceGUI application for configuring, running, and visualizing
1D oxygen diffusion-reaction simulations.
Supports multiple simulations with comparison visualization.
"""

from nicegui import ui, run, events
import asyncio
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from dataclasses import dataclass, field
from typing import Optional
import json
from datetime import datetime

from culturemods.rxd_fipy_1d import SimulationConfig, SimulationResult, run_simulation
from culturemods import kinetics

#FIXME - move to utility module
def calculate_discrete_rate_at_heights(df):
    time_col = 't_s'
    height_col = 'z_idx'
    df = df.sort_values([height_col, time_col])  # Important: sort first!

    # Calculate derivative each per-height time series
    dt = df.groupby(height_col)['t_s'].diff()
    dC = df.groupby(height_col)['C'].diff()
    df['discrete_rate'] = -60 * dC / dt #umolar per min decrease
    return df


# Color palette for multiple simulations
COLORS = [
    '#00d4aa',  # teal
    '#ff6b6b',  # red
    '#4ecdc4',  # cyan
    '#ffe66d',  # yellow
    '#95e1d3',  # mint
    '#f38181',  # coral
    '#aa96da',  # purple
    '#fcbad3',  # pink
]

# Line styles available
LINE_STYLES = {
    'solid': None,  # plotly default
    'dash': 'dash',
    'dot': 'dot',
    'dashdot': 'dashdot',
}


@dataclass
class SimulationEntry:
    """A single simulation with its configuration and results."""
    name: str
    config_values: dict
    result: Optional[SimulationResult] = None
    df: Optional[object] = None  # pandas DataFrame
    visible: bool = True
    color: str = '#00d4aa'
    line_style: str = 'solid'


class SimulationGUI:
    """GUI for oxygen diffusion simulation with multi-simulation support."""

    def __init__(self):
        # List of simulation entries
        self.simulations: list[SimulationEntry] = []
        self.next_sim_id = 1

        # Currently editing simulation index (-1 = adding new)
        self.editing_index: int = -1

        # Editor config values (for the form)
        self.editor_values = self._default_config()

        # Editor style values
        self.editor_color = '#00d4aa'
        self.editor_line_style = 'solid'

        # Visualization parameters
        self.probe_height_mm = 1.0
        self.profile_step = 0
        self.heatmap_sim_name = None  # Which simulation to use for heatmap (None = first visible)
        self.timeseries_mode = 'concentration'  # 'concentration' or 'rate'

        # Theme
        self.dark_mode = True

        # UI element references
        self.timeseries_plot = None
        self.profile_plot = None
        self.status_label = None
        self.probe_slider = None
        self.step_slider = None
        self.heatmap_select = None
        self.progress_bar = None
        self.sim_list_container = None
        self.editor_container = None
        self.editor_card = None
        self.dark_mode_toggle = None

        # Progress tracking
        self.current_step = 0
        self.total_steps = 0
        self.is_running = False

    def _default_config(self) -> dict:
        """Return default configuration values."""
        return {
            # Mode selection
            'mode': 'suspension',  # 'suspension' or 'monolayer'

            # Common parameters
            'D': 3e-3,
            'C_air': 200.0,

            # Suspension mode parameters
            'nz': 100,
            'rate': 10,  # umolar / min
            'k1': 0.0,
            'flux_bottom': 0.0,  # fmol/mm²/s
            'top_constraint': 'open',
            'dt': 1.0,
            'steps': 1800,
            'C_initial_fraction': 1.0,
            'halt_on_C_zero': True,

            # Monolayer mode parameters
            'flux': 100.0,  # fmol/mm²/s (OCR)
            'media_vol': 100.0,  # µL
            'well_radius': 3.2,  # mm (96-well default)
            'duration_hrs': 4.0,  # simulation duration in hours
            'time_step_s': 30,  # time step in seconds
        }

    def _get_next_color(self) -> str:
        """Get the next color from the palette."""
        used_colors = {s.color for s in self.simulations}
        for color in COLORS:
            if color not in used_colors:
                return color
        # If all colors used, cycle through
        return COLORS[len(self.simulations) % len(COLORS)]

    def build_ui(self):
        """Build the main UI layout."""
        self.dark = ui.dark_mode()
        self.dark.enable()

        with ui.header().classes('items-center justify-between'):
            ui.label('Culture Well O₂ Diffusion-Reaction Simulator').classes('text-2xl font-bold')
            with ui.row().classes('items-center gap-2'):
                ui.icon('light_mode').classes('text-yellow-400')
                self.dark_mode_toggle = ui.switch(
                    value=self.dark_mode,
                    on_change=self._toggle_theme
                )
                ui.icon('dark_mode').classes('text-blue-400')

        with ui.row().classes('w-full gap-4 p-4'):
            # Left panel - Simulations list and editor
            with ui.column().classes('w-96 gap-4'):
                self._build_simulations_panel()

            # Right panel - Visualization
            with ui.column().classes('flex-grow gap-4'):
                self._build_visualization_panel()

    def _build_simulations_panel(self):
        """Build the simulations list and editor panel."""
        # Simulations list card
        with ui.card().classes('w-full'):
            with ui.row().classes('items-center justify-between mb-2'):
                ui.label('Simulations').classes('text-lg font-bold')
                ui.button(icon='add', on_click=self._add_new_simulation).props('flat dense').tooltip('New simulation')
            with ui.row().classes('items-center justify-between mb-2'):
                ui.upload(
                    on_upload=self._load_simulation,
                    auto_upload=True,
                    multiple=True,
                    label = "Upload saved",
                ).props('flat dense accept=.json').classes('w-full').tooltip('Load simulations from file(s)')

            self.sim_list_container = ui.column().classes('w-full gap-1')
            self._refresh_sim_list()

            # Status and progress
            self.progress_bar = ui.linear_progress(value=0, show_value=False).classes('w-full mt-2')
            self.progress_bar.visible = False
            self.status_label = ui.label('Ready').classes('text-sm text-gray-400')

        # Editor card (initially hidden)
        self.editor_card = ui.card().classes('w-full')
        self.editor_card.visible = False
        with self.editor_card:
            self.editor_container = ui.column().classes('w-full gap-2')

    def _refresh_sim_list(self):
        """Refresh the simulation list display."""
        self.sim_list_container.clear()
        with self.sim_list_container:
            if not self.simulations:
                ui.label('No simulations yet. Click + to add one.').classes('text-sm text-gray-400 italic')
            else:
                for i, sim in enumerate(self.simulations):
                    self._build_sim_list_item(i, sim)

    def _build_sim_list_item(self, index: int, sim: SimulationEntry):
        """Build a single simulation list item."""
        with ui.row().classes('w-full items-center gap-2 p-2 rounded hover:bg-gray-800'):
            # Color indicator and visibility toggle
            ui.checkbox(
                value=sim.visible,
                on_change=lambda e, idx=index: self._toggle_visibility(idx, e.value)
            ).props('dense').style(f'color: {sim.color}')

            # Color dot
            ui.element('div').classes('w-3 h-3 rounded-full').style(f'background-color: {sim.color}')

            # Name and status
            with ui.column().classes('flex-grow gap-0'):
                ui.label(sim.name).classes('text-sm font-medium')
                if sim.result:
                    steps = len(sim.df['step'].unique()) if sim.df is not None else 0
                    ui.label(f'{steps} steps').classes('text-xs text-gray-400')
                else:
                    ui.label('Not run').classes('text-xs text-yellow-400')

            # Action buttons
            with ui.row().classes('gap-1'):
                ui.button(icon='edit', on_click=lambda idx=index: self._edit_simulation(idx)).props('flat dense size=sm').tooltip('Edit')
                ui.button(icon='play_arrow', on_click=lambda idx=index: self._run_single_simulation(idx)).props('flat dense size=sm').tooltip('Run')
                ui.button(icon='download', on_click=lambda idx=index: self._save_single_simulation(idx)).props('flat dense size=sm').tooltip('Save to JSON')
                ui.button(icon='table_chart', on_click=lambda idx=index: self._export_csv(idx)).props('flat dense size=sm').tooltip('Export CSV')
                ui.button(icon='delete', on_click=lambda idx=index: self._delete_simulation(idx)).props('flat dense size=sm color=red').tooltip('Delete')

    def _toggle_visibility(self, index: int, visible: bool):
        """Toggle simulation visibility."""
        self.simulations[index].visible = visible
        self._update_timeseries_plot()
        self._update_profile_plot()

    def _toggle_theme(self, e):
        """Toggle between dark and light theme."""
        self.dark_mode = e.value
        if self.dark_mode:
            self.dark.enable()
        else:
            self.dark.disable()
        # Refresh plots with new theme
        self._update_timeseries_plot()
        self._update_profile_plot()

    def _get_plot_template(self) -> str:
        """Get the plotly template based on current theme."""
        return 'plotly_dark' if self.dark_mode else 'plotly_white'

    def _add_new_simulation(self):
        """Start adding a new simulation."""
        self.editing_index = -1
        self.editor_values = self._default_config()
        self.editor_color = self._get_next_color()
        self.editor_line_style = 'solid'
        self._show_editor(f'Simulation {self.next_sim_id}')

    def _edit_simulation(self, index: int):
        """Edit an existing simulation."""
        self.editing_index = index
        sim = self.simulations[index]
        self.editor_values = sim.config_values.copy()
        self.editor_color = sim.color
        self.editor_line_style = sim.line_style
        self._show_editor(sim.name, is_edit=True)

    def _show_editor(self, name: str, is_edit: bool = False):
        """Show the configuration editor."""
        self.editor_card.visible = True
        self.editor_container.clear()

        # Ensure mode exists in editor_values (for backwards compatibility)
        if 'mode' not in self.editor_values:
            self.editor_values['mode'] = 'suspension'

        with self.editor_container:
            # Header
            with ui.row().classes('w-full items-center justify-between mb-2'):
                ui.label('Edit Simulation' if is_edit else 'New Simulation').classes('text-lg font-bold')
                ui.button(icon='close', on_click=self._hide_editor).props('flat dense')

            # Name input
            name_input = ui.input('Name', value=name).classes('w-full')

            ui.separator()

            # Appearance
            ui.label('Appearance').classes('font-semibold text-blue-400')

            with ui.row().classes('w-full items-center gap-4'):
                ui.label('Color:').classes('text-sm')
                ui.color_input(
                    value=self.editor_color,
                    on_change=lambda e: setattr(self, 'editor_color', e.value)
                ).classes('w-24')

                ui.label('Line style:').classes('text-sm')
                ui.select(
                    list(LINE_STYLES.keys()),
                    value=self.editor_line_style,
                    on_change=lambda e: setattr(self, 'editor_line_style', e.value)
                ).classes('w-28')

            ui.separator()

            # Mode selection
            ui.label('Simulation Mode').classes('font-semibold text-blue-400')

            # Container for mode-specific parameters (created later)
            params_container = None

            def update_mode(e):
                self._update_editor('mode', e.value)
                self._rebuild_mode_params(params_container, name_input)

            ui.toggle(
                ['suspension', 'monolayer'],
                value=self.editor_values.get('mode', 'suspension'),
                on_change=update_mode
            ).classes('w-full')

            # Now create the params container after the toggle
            params_container = ui.column().classes('w-full gap-2')

            # Build initial mode parameters
            self._rebuild_mode_params(params_container, name_input)

    def _rebuild_mode_params(self, container, name_input):
        """Rebuild the mode-specific parameters in the editor."""
        container.clear()
        mode = self.editor_values.get('mode', 'suspension')

        with container:
            ui.separator()

            # Common parameters
            ui.label('Physical Parameters').classes('font-semibold text-blue-400')

            ui.number(
                'Diffusion coeff D (mm²/s)',
                value=self.editor_values['D'],
                format='%.2e',
                step=1e-4,
                on_change=lambda e: self._update_editor('D', e.value)
            ).classes('w-full')

            ui.number(
                'C_air / C_initial (µM)',
                value=self.editor_values['C_air'],
                format='%.1f',
                on_change=lambda e: self._update_editor('C_air', e.value)
            ).classes('w-full')

            if mode == 'suspension':
                self._build_suspension_params()
            else:
                self._build_monolayer_params()

            ui.separator()

            # Action buttons
            with ui.row().classes('w-full gap-2'):
                ui.button(
                    'Save & Run',
                    on_click=lambda: self._save_and_run(name_input.value)
                ).classes('flex-grow').props('color=primary')
                ui.button(
                    'Save',
                    on_click=lambda: self._save_simulation(name_input.value)
                ).classes('flex-grow')
                ui.button('Cancel', on_click=self._hide_editor).classes('flex-grow')

    def _build_suspension_params(self):
        """Build suspension mode parameters."""
        ui.number(
            'Media volume (µL)',
            value=self.editor_values.get('media_vol', 100.0),
            format='%.1f',
            step=10, min=10,
            on_change=lambda e: self._update_editor('media_vol', e.value)
        ).classes('w-full')

        ui.number(
            'Well radius (mm)',
            value=self.editor_values.get('well_radius', 3.2),
            format='%.2f',
            step=0.1, min=0.5,
            on_change=lambda e: self._update_editor('well_radius', e.value)
        ).classes('w-full')

        # Show calculated media height
        vol = self.editor_values.get('media_vol', 100.0)
        radius = self.editor_values.get('well_radius', 3.2)
        height = kinetics.media_vol_to_height(vol, radius)
        ui.label(f'Media height: {height:.2f} mm').classes('text-sm text-gray-400')

        ui.number(
            'Mesh points (nz)',
            value=self.editor_values['nz'],
            format='%.0f',
            min=10, max=500, step=10,
            on_change=lambda e: self._update_editor('nz', int(e.value))
        ).classes('w-full')

        ui.separator()

        ui.label('Reaction Parameters').classes('font-semibold text-blue-400')

        ui.number(
            'Zero-order rate (µM/min)',
            value=self.editor_values['rate'],
            format='%.1f',
            step=0.1,
            on_change=lambda e: self._update_editor('rate', e.value)
        ).classes('w-full')

        ui.number(
            'First-order rate k1',
            value=self.editor_values['k1'],
            format='%.6f',
            step=0.001, min=0,
            on_change=lambda e: self._update_editor('k1', e.value)
        ).classes('w-full')

        ui.number(
            'Bottom flux (fmol/mm²/s)',
            value=self.editor_values['flux_bottom'],
            format='%.2f',
            step=0.1, min=0,
            on_change=lambda e: self._update_editor('flux_bottom', e.value)
        ).classes('w-full')

        ui.select(
            ['open', 'sealed'],
            value=self.editor_values['top_constraint'],
            label='Top boundary',
            on_change=lambda e: self._update_editor('top_constraint', e.value)
        ).classes('w-full')

        ui.separator()

        ui.label('Time Parameters').classes('font-semibold text-blue-400')

        ui.number(
            'Time step dt (s)',
            value=self.editor_values['dt'],
            format='%.2f',
            step=0.1, min=0.01,
            on_change=lambda e: self._update_editor('dt', e.value)
        ).classes('w-full')

        ui.number(
            'Number of steps',
            value=self.editor_values['steps'],
            format='%.0f',
            min=1, max=10000, step=10,
            on_change=lambda e: self._update_editor('steps', int(e.value))
        ).classes('w-full')

        ui.number(
            'Initial C fraction',
            value=self.editor_values['C_initial_fraction'],
            format='%.2f',
            min=0, max=1, step=0.1,
            on_change=lambda e: self._update_editor('C_initial_fraction', e.value)
        ).classes('w-full')

        ui.checkbox(
            'Halt when C=0',
            value=self.editor_values['halt_on_C_zero'],
            on_change=lambda e: self._update_editor('halt_on_C_zero', e.value)
        )

    def _build_monolayer_params(self):
        """Build monolayer mode parameters."""
        ui.separator()

        ui.label('Monolayer Parameters').classes('font-semibold text-blue-400')

        ui.number(
            'O₂ flux / OCR (fmol/mm²/s)',
            value=self.editor_values.get('flux', 100.0),
            format='%.1f',
            step=10, min=0,
            on_change=lambda e: self._update_editor('flux', e.value)
        ).classes('w-full')

        ui.number(
            'Media volume (µL)',
            value=self.editor_values.get('media_vol', 100.0),
            format='%.1f',
            step=10, min=10,
            on_change=lambda e: self._update_editor('media_vol', e.value)
        ).classes('w-full')

        ui.number(
            'Well radius (mm)',
            value=self.editor_values.get('well_radius', 3.2),
            format='%.2f',
            step=0.1, min=0.5,
            on_change=lambda e: self._update_editor('well_radius', e.value)
        ).classes('w-full')

        # Show calculated media height
        vol = self.editor_values.get('media_vol', 100.0)
        radius = self.editor_values.get('well_radius', 3.2)
        height = kinetics.media_vol_to_height(vol, radius)
        ui.label(f'Media height: {height:.2f} mm').classes('text-sm text-gray-400')

        ui.separator()

        ui.label('Time Parameters').classes('font-semibold text-blue-400')

        ui.number(
            'Duration (hours)',
            value=self.editor_values.get('duration_hrs', 4.0),
            format='%.1f',
            step=0.5, min=0.1,
            on_change=lambda e: self._update_editor('duration_hrs', e.value)
        ).classes('w-full')

        ui.number(
            'Time step (s)',
            value=self.editor_values.get('time_step_s', 30),
            format='%.0f',
            step=10, min=1,
            on_change=lambda e: self._update_editor('time_step_s', int(e.value))
        ).classes('w-full')

    def _hide_editor(self):
        """Hide the configuration editor."""
        self.editor_card.visible = False
        self.editing_index = -1

    def _update_editor(self, key: str, value):
        """Update editor configuration value."""
        self.editor_values[key] = value

    def _save_simulation(self, name: str):
        """Save the current editor configuration."""
        if self.editing_index >= 0:
            # Update existing
            sim = self.simulations[self.editing_index]

            # Check if config values changed (not just appearance)
            config_changed = sim.config_values != self.editor_values

            # Update appearance (always safe)
            sim.name = name
            sim.color = self.editor_color
            sim.line_style = self.editor_line_style

            # Only clear results if config actually changed
            if config_changed:
                sim.config_values = self.editor_values.copy()
                sim.result = None
                sim.df = None
        else:
            # Create new
            sim = SimulationEntry(
                name=name,
                config_values=self.editor_values.copy(),
                color=self.editor_color,
                line_style=self.editor_line_style
            )
            self.simulations.append(sim)
            self.next_sim_id += 1

        self._hide_editor()
        self._refresh_sim_list()
        self._update_timeseries_plot()
        self._update_profile_plot()

    async def _save_and_run(self, name: str):
        """Save and immediately run the simulation."""
        self._save_simulation(name)
        # Run the just-saved simulation
        index = self.editing_index if self.editing_index >= 0 else len(self.simulations) - 1
        await self._run_single_simulation(index)

    def _delete_simulation(self, index: int):
        """Delete a simulation."""
        del self.simulations[index]
        self._refresh_sim_list()
        self._update_heatmap_select_options()
        self._update_timeseries_plot()
        self._update_profile_plot()

    def _save_single_simulation(self, index: int):
        """Save a single simulation to JSON file."""
        sim = self.simulations[index]
        data = {
            'name': sim.name,
            'config_values': sim.config_values,
            'color': sim.color,
            'line_style': sim.line_style,
            'result': sim.result.__getstate__() if sim.result else None,
        }
        json_str = json.dumps(data, indent=2)

        # Sanitize name for filename
        safe_name = ''.join(c if c.isalnum() or c in '-_' else '_' for c in sim.name)
        filename = f'{safe_name}.json'
        ui.download(json_str.encode('utf-8'), filename)

        self.status_label.text = f'Saved: {sim.name}'
        self.status_label.classes('text-green-400', remove='text-gray-400 text-yellow-400')

    def _export_csv(self, index: int):
        """Export simulation results as CSV file."""
        sim = self.simulations[index]
        if sim.df is None:
            self.status_label.text = f'No results to export for {sim.name}'
            self.status_label.classes('text-yellow-400', remove='text-gray-400 text-green-400 text-red-400')
            return

        csv_str = sim.df.to_csv(index=False)

        # Sanitize name for filename
        safe_name = ''.join(c if c.isalnum() or c in '-_' else '_' for c in sim.name)
        filename = f'{safe_name}.csv'
        ui.download(csv_str.encode('utf-8'), filename)

        self.status_label.text = f'Exported CSV: {sim.name}'
        self.status_label.classes('text-green-400', remove='text-gray-400 text-yellow-400')

    async def _load_simulation(self, e: events.UploadEventArguments):
        """Load a simulation from JSON file."""
        try:
            content = await e.file.text()
            data = json.loads(content)

            # Create simulation entry
            sim = SimulationEntry(
                name=data['name'],
                config_values=data['config_values'],
                color=data.get('color', self._get_next_color()),
                line_style=data.get('line_style', 'solid'),
            )

            # Restore result if present
            if data.get('result'):
                sim.result = SimulationResult.__new__(SimulationResult)
                sim.result.__setstate__(data['result'])
                sim.df = sim.result.to_dataframe()
                #FIXME - support for old dataframes without discrete_rate
                #FIXME - this shouldn't be done in gui module
                sim.df = calculate_discrete_rate_at_heights(sim.df)

            self.simulations.append(sim)
            self.next_sim_id += 1

            self._refresh_sim_list()
            self._update_slider_ranges()
            self._update_timeseries_plot()
            self._update_profile_plot()

            has_result = 'with results' if sim.result else 'no results'
            self.status_label.text = f'Loaded: {sim.name} ({has_result})'
            self.status_label.classes('text-green-400', remove='text-gray-400 text-yellow-400 text-red-400')

        except Exception as ex:
            self.status_label.text = f'Load error: {str(ex)}'
            self.status_label.classes('text-red-400', remove='text-gray-400 text-yellow-400 text-green-400')
            raise

    def _create_config_from_values(self, values: dict) -> SimulationConfig:
        """Create SimulationConfig from config values dict."""
        C_air = values['C_air']

        # Calculate media height from volume and well radius
        media_vol = values.get('media_vol', 100.0)
        well_radius = values.get('well_radius', 3.2)
        L = kinetics.media_vol_to_height(media_vol, well_radius)

        # Convert volumetric reaction rate - umolar/min to umolar/s
        rate_umolar_per_s = values['rate']/60
        k = rate_umolar_per_s / C_air

        dt = values['dt']
        # Convert bottom flux
        flux_bottom = values['flux_bottom'] / C_air * dt

        return SimulationConfig(
            D=values['D'],
            C_air=C_air,
            L=L,
            nz=values['nz'],
            k=k,
            k1=values['k1'],
            flux_bottom=flux_bottom,
            top_constraint=values['top_constraint'],
            dt=dt,
            steps=values['steps'],
            C_initial_fraction=values['C_initial_fraction'],
            halt_on_C_zero=values['halt_on_C_zero'],
        )

    async def _run_single_simulation(self, index: int):
        """Run a single simulation by index."""
        if self.is_running:
            return

        sim = self.simulations[index]
        mode = sim.config_values.get('mode', 'suspension')

        self.is_running = True
        self.status_label.text = f'Running {sim.name}...'
        self.status_label.classes('text-yellow-400', remove='text-gray-400 text-green-400 text-red-400')

        self.progress_bar.visible = True
        self.progress_bar.value = 0
        self.current_step = 0

        try:
            if mode == 'monolayer':
                await self._run_monolayer_simulation(sim)
            else:
                await self._run_suspension_simulation(sim)

            # Update slider ranges based on all simulations
            self._update_slider_ranges()

            self._refresh_sim_list()
            self._update_timeseries_plot()
            self._update_profile_plot()

            actual_steps = len(sim.df['step'].unique()) if sim.df is not None else 0
            self.status_label.text = f'{sim.name}: {actual_steps} steps'
            self.status_label.classes('text-green-400', remove='text-gray-400 text-yellow-400 text-red-400')

        except Exception as e:
            self.status_label.text = f'Error: {str(e)}'
            self.status_label.classes('text-red-400', remove='text-gray-400 text-yellow-400 text-green-400')

        finally:
            self.is_running = False
            await asyncio.sleep(1.0)
            self.progress_bar.visible = False

    async def _run_suspension_simulation(self, sim: SimulationEntry):
        """Run a suspension mode simulation using FiPy."""
        config = self._create_config_from_values(sim.config_values)
        self.total_steps = config.steps

        progress_timer = ui.timer(0.1, self._update_progress)

        sim.result = await run.io_bound(
            run_simulation,
            config,
            record_every=1,
            verbose=False,
            step_callback=self._step_callback
        )

        progress_timer.cancel()
        self.progress_bar.value = 1.0

        sim.df = sim.result.to_dataframe()
        #FIXME - this shouldn't be done in gui module
        sim.df = calculate_discrete_rate_at_heights(sim.df)

    async def _run_monolayer_simulation(self, sim: SimulationEntry):
        """Run a monolayer mode simulation using kinetics.py analytical solution."""
        import pandas as pd

        values = sim.config_values
        flux_fmol = values.get('flux', 100.0)
        media_vol = values.get('media_vol', 100.0)
        well_radius = values.get('well_radius', 3.2)
        c_initial = values.get('C_air', 200.0)
        duration_hrs = values.get('duration_hrs', 4.0)
        time_step_s = values.get('time_step_s', 30)
        D = values.get('D', 3e-3)

        # Convert flux to units expected by kinetics.py (µmol/mm²/s)
        Q = kinetics.flux_units_convert(flux_fmol)

        # Calculate media height from volume
        media_height = kinetics.media_vol_to_height(media_vol, well_radius)

        # Generate time points
        duration_s = int(duration_hrs * 3600)
        time_points = list(range(0, duration_s + 1, time_step_s))
        self.total_steps = len(time_points)

        # Number of spatial points
        nz = 20
        dz = media_height / nz

        progress_timer = ui.timer(0.1, self._update_progress)


        # Run simulation in background
        def run_kinetics():
            points = []
            for step, t_s in enumerate(time_points):
                self.current_step = step
                for z_idx in range(nz):
                    height = z_idx * dz
                    # kinetics.concentration expects position from bottom
                    c = kinetics.concentration(height, t_s, Q, c_initial, media_height)
                    c = max(0, c)  # Clamp to non-negative
                    c_star = c / c_initial

                    points.append({
                        'k': flux_fmol,  # Store flux as 'k' for compatibility
                        'step': step,
                        'c_star': c_star,
                        'z_idx': z_idx,
                        'C': c,
                        'height_mm': height,
                        't_hrs': t_s / 3600,
                        't_mins': t_s / 60,
                        't_s': t_s,
                    })
            return points

        points = await run.io_bound(run_kinetics)

        progress_timer.cancel()
        self.progress_bar.value = 1.0

        # Create DataFrame directly (no SimulationResult for monolayer)
        sim.df = pd.DataFrame(points)

        #FIXME - this shouldn't be done in gui module
        sim.df = calculate_discrete_rate_at_heights(sim.df)

        #FIXME - unify this simulation result with SimulationResult class fromr xd_fipy_1d in
        # a cleaner way

        # Create a minimal result object for compatibility
        # Store config info needed for plotting
        sim.result = type('MonolayerResult', (), {
            'config': type('Config', (), {
                'L': media_height,
                'nz': nz,
                'dz': dz,
                'C_air': c_initial,
                'dt': time_step_s,
            })(),
            'final_profile': None,
            'points': points,
            'to_dataframe': lambda: sim.df,
            'get_profile_at_step': lambda _, step: sim.df[sim.df['step'] == step].sort_values('z_idx'),
            '__getstate__': lambda: {
                'config': {
                    'mode': 'monolayer',
                    'flux': flux_fmol,
                    'media_vol': media_vol,
                    'well_radius': well_radius,
                    'C_air': c_initial,
                    'D': D,
                    'duration_hrs': duration_hrs,
                    'time_step_s': time_step_s,
                },
                'points': points,
                'final_profile': None,
            },
        })()

    def _step_callback(self, step: int):
        """Callback for simulation progress."""
        self.current_step = step

    def _update_progress(self):
        """Update progress bar."""
        if self.total_steps > 0:
            self.progress_bar.value = self.current_step / self.total_steps

    def _update_slider_ranges(self):
        """Update slider ranges based on all simulations."""
        max_L = 3.1
        max_steps = 100

        for sim in self.simulations:
            if sim.result:
                max_L = max(max_L, sim.result.config.L)
                if sim.df is not None:
                    max_steps = max(max_steps, sim.df['step'].max())

        if self.probe_slider:
            self.probe_slider._props['max'] = max_L
        if self.step_slider:
            self.step_slider._props['max'] = max_steps

        # Update heatmap simulation selector
        self._update_heatmap_select_options()

    def _build_visualization_panel(self):
        """Build the visualization panel with plots."""
        # Timeseries plot card
        with ui.card().classes('w-full'):
            with ui.row().classes('items-center gap-4 mb-2'):
                self.timeseries_title = ui.label('O₂ @ Probe Height vs Time').classes('text-lg font-bold')
                ui.toggle(
                    {'concentration': 'Conc', 'rate': 'dC/dt'},
                    value=self.timeseries_mode,
                    on_change=self._on_timeseries_mode_change
                ).props('dense')
                ui.label('Probe height:').classes('text-sm')
                self.probe_slider = ui.slider(
                    min=0, max=3.1, step=0.1, value=self.probe_height_mm,
                    on_change=self._on_probe_height_change
                ).classes('w-48')
                self.probe_label = ui.label(f'{self.probe_height_mm:.1f} mm').classes('text-sm w-16')

            self.timeseries_plot = ui.plotly({}).classes('w-full h-[400px]')

        # Profile plot card
        with ui.card().classes('w-full'):
            with ui.row().classes('items-center gap-4 mb-2'):
                ui.label('Vertical Concentration Profile').classes('text-lg font-bold')
                ui.label('Time step:').classes('text-sm')
                self.step_slider = ui.slider(
                    min=0, max=1800, step=60, value=self.profile_step,
                    on_change=self._on_step_change
                ).classes('w-48')
                self.step_label = ui.label('Step 0 (0.0 s)').classes('text-sm w-32')
            with ui.row().classes('items-center gap-4 mb-2'):
                ui.label('Gradient Profile:').classes('text-sm')
                self.heatmap_select = ui.select(
                    options=[],
                    value=None,
                    on_change=self._on_heatmap_sim_change
                ).classes('w-40').props('dense')

            self.profile_plot = ui.plotly({}).classes('w-full h-[400px]')

        self._update_timeseries_plot()
        self._update_profile_plot()

    def _on_probe_height_change(self, e):
        """Handle probe height change."""
        self.probe_height_mm = e.value
        self.probe_label.text = f'{self.probe_height_mm:.1f} mm'
        self._update_timeseries_plot()

    def _on_timeseries_mode_change(self, e):
        """Handle timeseries mode toggle between concentration and rate."""
        self.timeseries_mode = e.value
        # Update title
        if self.timeseries_mode == 'concentration':
            self.timeseries_title.text = 'O₂ @ Probe Height vs Time'
        else:
            self.timeseries_title.text = 'O₂ Change @ Probe Height vs Time'
        self._update_timeseries_plot()

    def _on_step_change(self, e):
        """Handle step slider change."""
        self.profile_step = int(e.value)
        # Update label with time from first visible simulation
        for sim in self.simulations:
            if sim.visible and sim.result:
                t_s = self.profile_step * sim.result.config.dt
                self.step_label.text = f'Step {self.profile_step} ({t_s:.1f} s)'
                break
        self._update_profile_plot()

    def _on_heatmap_sim_change(self, e):
        """Handle heatmap simulation selection change."""
        self.heatmap_sim_name = e.value
        self._update_profile_plot()

    def _update_heatmap_select_options(self):
        """Update the heatmap simulation selector options."""
        if not self.heatmap_select:
            return

        # Get list of simulations with results
        options = {sim.name: sim.name for sim in self.simulations if sim.result is not None}

        self.heatmap_select.options = options
        self.heatmap_select.update()

        # If current selection is no longer valid, reset to None (first visible)
        if self.heatmap_sim_name and self.heatmap_sim_name not in options:
            self.heatmap_sim_name = None
            self.heatmap_select.value = None

    def _update_timeseries_plot(self):
        """Update the timeseries plot with all visible simulations."""
        fig = go.Figure()
        template = self._get_plot_template()
        show_rate = self.timeseries_mode == 'rate'

        visible_sims = [s for s in self.simulations if s.visible and s.df is not None]

        y_label = 'dC/dt (µM/min)' if show_rate else 'O₂ Concentration (µM)'
        y_col = 'discrete_rate' if show_rate else 'C'

        if not visible_sims:
            fig.update_layout(
                template=template,
                margin=dict(l=60, r=20, t=30, b=50),
                xaxis_title='Time (minutes)',
                yaxis_title=y_label,
                yaxis=dict(range=[0, 220] if not show_rate else None),
            )
            fig.add_annotation(
                text="Add and run simulations to see results",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=14, color="gray")
            )
            self.timeseries_plot.update_figure(fig)
            return

        y_max = 0
        y_min = 0
        c_air_max = 0

        for sim in visible_sims:
            config = sim.result.config
            probe_idx = int(self.probe_height_mm / config.dz)
            probe_idx = max(0, min(probe_idx, config.nz - 1))

            df_probe = sim.df[sim.df['z_idx'] == probe_idx].sort_values('t_mins')

            # Get line dash style
            dash = LINE_STYLES.get(sim.line_style)

            # Get y values, handling missing discrete_rate column
            if y_col in df_probe.columns:
                y_values = df_probe[y_col]
            else:
                y_values = df_probe['C']

            fig.add_trace(go.Scatter(
                x=df_probe['t_mins'],
                y=y_values,
                mode='lines',
                name=sim.name,
                line=dict(color=sim.color, width=2, dash=dash)
            ))

            if y_col in sim.df.columns:
                col_values = sim.df[y_col].dropna()
                if len(col_values) > 0:
                    y_max = max(y_max, col_values.max())
                    y_min = min(y_min, col_values.min())
            c_air_max = max(c_air_max, config.C_air)

        # Set y-axis range
        if show_rate:
            # For rate, center around zero with some padding
            y_range = [y_min * 1.1 if y_min < 0 else y_min - 0.1,
                       y_max * 1.1 if y_max > 0 else 0.1]
        else:
            y_range = [-0.1, max(y_max, c_air_max) * 1.1]

        fig.update_layout(
            template=template,
            margin=dict(l=60, r=20, t=30, b=50),
            xaxis_title='Time (minutes)',
            yaxis_title=y_label,
            yaxis=dict(range=y_range),
            showlegend=True,
            legend=dict(x=0.02, y=0.98),
        )

        # Add C_air reference line only for concentration mode
        if not show_rate:
            fig.add_hline(y=c_air_max, line_dash="dash", line_color="gray",
                          annotation_text=f"C_air = {c_air_max} µM")

        self.timeseries_plot.update_figure(fig)

    def _update_profile_plot(self):
        """Update the vertical profile plot with all visible simulations."""
        template = self._get_plot_template()

        fig = make_subplots(
            rows=1, cols=2,
            column_widths=[0.15, 0.85],
            horizontal_spacing=0.02,
            shared_yaxes=True
        )

        visible_sims = [s for s in self.simulations if s.visible and s.result is not None]

        if not visible_sims:
            fig.update_layout(
                template=template,
                margin=dict(l=60, r=20, t=30, b=50),
            )
            fig.add_annotation(
                text="Add and run simulations to see results",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=14, color="gray")
            )
            self.profile_plot.update_figure(fig)
            return

        c_max = 0
        max_L = 0

        # Find simulation to use for heatmap
        heatmap_sim = None
        if self.heatmap_sim_name:
            # Use selected simulation if it exists and has results
            for sim in self.simulations:
                if sim.name == self.heatmap_sim_name and sim.result is not None:
                    heatmap_sim = sim
                    break
        if heatmap_sim is None:
            # Default to first visible simulation
            heatmap_sim = visible_sims[0]

        heatmap_config = heatmap_sim.result.config
        heatmap_profile = heatmap_sim.result.get_profile_at_step(self.profile_step)

        if len(heatmap_profile) > 0:
            heights = heatmap_profile['height_mm'].values
            concentrations = heatmap_profile['C'].values
            z_data = np.array(concentrations).reshape(-1, 1)

            fig.add_trace(go.Heatmap(
                z=z_data,
                x=[0],
                y=heights,
                colorscale='RdYlBu_r',
                zmin=0,
                zmax=heatmap_config.C_air,
                showscale=False,
                hovertemplate=f'{heatmap_sim.name}<br>Height: %{{y:.2f}} mm<br>C: %{{z:.1f}} µM<extra></extra>'
            ), row=1, col=1)

        # Add line plots for all visible simulations
        for sim in visible_sims:
            config = sim.result.config
            profile_df = sim.result.get_profile_at_step(self.profile_step)

            if len(profile_df) == 0:
                continue

            heights = profile_df['height_mm'].values
            concentrations = profile_df['C'].values

            # Get line dash style
            dash = LINE_STYLES.get(sim.line_style)

            fig.add_trace(go.Scatter(
                x=concentrations,
                y=heights,
                mode='lines',
                name=sim.name,
                line=dict(color=sim.color, width=2, dash=dash),
                hovertemplate=f'{sim.name}<br>C: %{{x:.1f}} µM<br>Height: %{{y:.2f}} mm<extra></extra>'
            ), row=1, col=2)

            c_max = max(c_max, np.max(concentrations))
            max_L = max(max_L, config.L)

            # Add probe marker
            probe_idx = int(self.probe_height_mm / config.dz)
            probe_idx = max(0, min(probe_idx, len(concentrations) - 1))
            if probe_idx < len(concentrations):
                fig.add_trace(go.Scatter(
                    x=[concentrations[probe_idx]],
                    y=[heights[probe_idx]],
                    mode='markers',
                    marker=dict(color=sim.color, size=8, symbol='diamond'),
                    showlegend=False,
                    hovertemplate=f'{sim.name} Probe<br>C: %{{x:.1f}} µM<extra></extra>'
                ), row=1, col=2)

        fig.update_layout(
            template=template,
            margin=dict(l=60, r=20, t=30, b=50),
            showlegend=True,
            legend=dict(x=1.02, y=0.98),
        )

        fig.update_xaxes(showticklabels=False, row=1, col=1)
        fig.update_xaxes(title_text='O₂ Concentration (µM)', range=[0, c_max * 1.1], row=1, col=2)
        fig.update_yaxes(title_text='Height (mm)', row=1, col=1)

        # Annotations
        fig.add_annotation(
            x=0, y=max_L, xref='x', yref='y',
            text="Air", showarrow=False,
            font=dict(size=10), xanchor='center', yanchor='bottom'
        )
        fig.add_annotation(
            x=0, y=0, xref='x', yref='y',
            text="Bottom", showarrow=False,
            font=dict(size=10), xanchor='center', yanchor='top'
        )

        self.profile_plot.update_figure(fig)


@ui.page('/')
def index():
    """Main page."""
    app = SimulationGUI()
    app.build_ui()


def main():
    """Main entry point."""
    ui.run(title='O₂ Diffusion Simulator', port=8081, reload=False)


if __name__ in {'__main__', '__mp_main__'}:
    main()
