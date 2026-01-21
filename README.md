# culturemods

Oxygen diffusion-reaction models for cell culture systems.

## Installation

```bash
pip install -e .
```

For GUI support:
```bash
pip install -e ".[gui]"
```

For development:
```bash
pip install -e ".[dev,gui]"
```

## Usage

### CLI

```bash
# Plot O2 at bottom over time for different OCR values
culturemods o2-at-bottom-by-ocr 100 150 200

# Plot O2 at bottom over time for different media volumes
culturemods o2-at-bottom-by-vol 50 100 150
```

### GUI

```bash
culturemods-gui
```

### Python API

```python
from culturemods import concentration, flux_units_convert, media_vol_to_height
from culturemods import SimulationConfig, run_simulation

# Analytical kinetics (monolayer)
flux = flux_units_convert(100)  # fmol/mm²/s to µmol/mm²/s
height = media_vol_to_height(100)  # µL to mm
c = concentration(0, 3600, flux, c_initial=200, media_height=height)

# FiPy simulation (suspension)
config = SimulationConfig(k=0.5, steps=100)
result = run_simulation(config)
df = result.to_dataframe()
```
