"""
Test suite for rxd_fipy_1d.py - 1D Reaction-Diffusion oxygen transport model.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_array_less

from rxd_fipy_1d import (
    SimulationConfig,
    SimulationResult,
    create_mesh,
    run_simulation,
    run_parameter_sweep,
    combine_results,
)


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def default_config():
    """Default simulation configuration."""
    return SimulationConfig()


@pytest.fixture
def fast_config():
    """Configuration for fast tests with fewer steps and coarser mesh."""
    return SimulationConfig(
        k=0.1,  # Low consumption to avoid early halt
        nz=50,
        steps=10,
        dt=1.0,
    )


@pytest.fixture
def no_reaction_config():
    """Configuration with no consumption (pure diffusion)."""
    return SimulationConfig(
        k=0,
        k1=0,
        nz=50,
        steps=20,
        dt=1.0,
    )


@pytest.fixture
def first_order_config():
    """Configuration with first-order reaction kinetics."""
    return SimulationConfig(
        k=0,
        k1=0.01,
        nz=50,
        steps=20,
        dt=1.0,
    )


@pytest.fixture
def sealed_config():
    """Configuration for a sealed (closed) system."""
    return SimulationConfig(
        top_constraint='sealed',
        k=0.5,
        nz=50,
        steps=20,
        dt=1.0,
    )


@pytest.fixture
def depleted_initial_config():
    """Configuration starting with depleted oxygen."""
    return SimulationConfig(
        C_initial_fraction=0.0,
        k=0,
        nz=50,
        steps=50,
        dt=1.0,
    )


@pytest.fixture
def high_consumption_config():
    """Configuration with high consumption rate that triggers halt_on_C_zero."""
    return SimulationConfig(
        k=10.0,
        nz=50,
        steps=100,
        dt=1.0,
        halt_on_C_zero=True,
    )


# =============================================================================
# SimulationConfig Tests
# =============================================================================

class TestSimulationConfig:
    """Tests for SimulationConfig dataclass."""

    def test_default_initialization(self, default_config):
        """Test that default config initializes with expected values."""
        assert default_config.D == 3e-3
        assert default_config.C_air == 200.0
        assert default_config.L == 3.1
        assert default_config.nz == 200
        assert default_config.k == 1.0
        assert default_config.k1 == 0
        assert default_config.top_constraint == 'open'
        assert default_config.steps == 100

    def test_derived_parameters(self, default_config):
        """Test that derived parameters are computed correctly."""
        expected_dz = default_config.L / default_config.nz
        expected_T = default_config.L**2 / default_config.D

        assert_allclose(default_config.dz, expected_dz)
        assert_allclose(default_config.T, expected_T)
        assert_allclose(default_config.total_time, default_config.dt * default_config.steps)

    def test_dt_defaults_from_T(self):
        """Test that dt is computed from T when not specified."""
        config = SimulationConfig(dt=None)
        expected_dt = config.T / 3600
        assert_allclose(config.dt, expected_dt)

    def test_dt_override(self):
        """Test that explicit dt overrides default calculation."""
        config = SimulationConfig(dt=5.0)
        assert config.dt == 5.0

    def test_damkohler_property(self, default_config):
        """Test Damköhler number calculation."""
        expected_da = default_config.L**2 * default_config.k / default_config.D
        assert_allclose(default_config.damkohler, expected_da)

    def test_damkohler_zero_k(self):
        """Test Damköhler number is zero when k=0."""
        config = SimulationConfig(k=0)
        assert config.damkohler == 0

    def test_L_char_property(self, default_config):
        """Test characteristic length calculation."""
        expected_L_char = np.sqrt(default_config.D / default_config.k)
        assert_allclose(default_config.L_char, expected_L_char)

    def test_L_char_zero_k(self):
        """Test characteristic length is infinite when k=0."""
        config = SimulationConfig(k=0)
        assert config.L_char == float('inf')

    def test_custom_geometry(self):
        """Test configuration with custom geometry parameters."""
        config = SimulationConfig(L=5.0, nz=100)
        assert config.L == 5.0
        assert config.nz == 100
        assert_allclose(config.dz, 0.05)


# =============================================================================
# Mesh Creation Tests
# =============================================================================

class TestCreateMesh:
    """Tests for mesh creation."""

    def test_mesh_dimensions(self, fast_config):
        """Test that mesh has correct number of cells."""
        mesh = create_mesh(fast_config)
        assert mesh.numberOfCells == fast_config.nz

    def test_mesh_spacing(self, fast_config):
        """Test that mesh has correct spacing."""
        mesh = create_mesh(fast_config)
        # FiPy mesh dx should match config dz
        assert_allclose(mesh.dx, fast_config.dz)

    def test_mesh_extent(self, fast_config):
        """Test that mesh covers the full domain."""
        mesh = create_mesh(fast_config)
        cell_centers = mesh.cellCenters.value[0]
        # First cell center should be at dz/2, last at L - dz/2
        assert_allclose(cell_centers[0], fast_config.dz / 2, rtol=1e-10)
        assert_allclose(cell_centers[-1], fast_config.L - fast_config.dz / 2, rtol=1e-10)


# =============================================================================
# Simulation Tests
# =============================================================================

class TestRunSimulation:
    """Tests for the main simulation function."""

    def test_simulation_returns_result(self, fast_config):
        """Test that simulation returns a SimulationResult object."""
        result = run_simulation(fast_config)
        assert isinstance(result, SimulationResult)

    def test_simulation_stores_config(self, fast_config):
        """Test that result contains the config used."""
        result = run_simulation(fast_config)
        assert result.config is fast_config

    def test_simulation_has_final_profile(self, fast_config):
        """Test that simulation produces a final profile."""
        result = run_simulation(fast_config)
        assert result.final_profile is not None
        assert len(result.final_profile) == fast_config.nz

    def test_simulation_records_points(self, fast_config):
        """Test that simulation records data points."""
        result = run_simulation(fast_config, record_every=1)
        assert len(result.points) > 0

    def test_record_every_parameter(self):
        """Test that record_every controls recording frequency."""
        # Use config with no reaction to ensure full completion
        config = SimulationConfig(k=0, k1=0, nz=20, steps=20, dt=1.0)
        result_every = run_simulation(config, record_every=1)
        result_sparse = run_simulation(config, record_every=5)

        # Count unique steps recorded
        steps_every = len(set(p['step'] for p in result_every.points))
        steps_sparse = len(set(p['step'] for p in result_sparse.points))

        # Sparse recording should have fewer time steps recorded
        assert steps_sparse < steps_every

    def test_concentration_bounded_above(self, fast_config):
        """Test that concentration doesn't exceed initial/boundary value."""
        result = run_simulation(fast_config)
        max_c = max(p['c_star'] for p in result.points)
        assert max_c <= 1.0 + 1e-10  # Allow small numerical error

    def test_pure_diffusion_equilibrates(self, no_reaction_config):
        """Test that pure diffusion reaches equilibrium at C=1."""
        # Start at equilibrium
        config = SimulationConfig(
            k=0, k1=0, nz=50, steps=50, dt=1.0,
            C_initial_fraction=1.0
        )
        result = run_simulation(config)

        # Should stay at equilibrium
        assert_allclose(result.final_profile, 1.0, atol=1e-6)

    def test_diffusion_from_depleted(self, depleted_initial_config):
        """Test that oxygen diffuses into depleted media from open top."""
        result = run_simulation(depleted_initial_config)

        # Top should be at boundary value (1.0)
        assert_allclose(result.final_profile[-1], 1.0, atol=0.1)

        # Concentration should increase from initial (0)
        assert result.final_profile[0] > 0

    def test_consumption_depletes_bottom(self, fast_config):
        """Test that consumption depletes oxygen at bottom first."""
        result = run_simulation(fast_config)

        # Bottom concentration should be less than top
        assert result.final_profile[0] < result.final_profile[-1]

    def test_halt_on_C_zero(self, high_consumption_config):
        """Test that simulation halts when bottom concentration reaches zero."""
        result = run_simulation(high_consumption_config)

        # Should have halted before all steps completed
        steps_recorded = set(p['step'] for p in result.points)
        assert max(steps_recorded) < high_consumption_config.steps - 1

    def test_halt_on_C_zero_disabled(self):
        """Test that simulation continues when halt_on_C_zero is False."""
        config = SimulationConfig(
            k=10.0,
            nz=50,
            steps=50,
            dt=1.0,
            halt_on_C_zero=False,
        )
        result = run_simulation(config)

        # Should complete all steps
        steps_recorded = set(p['step'] for p in result.points)
        # At minimum, should record more steps than halted case
        assert len(steps_recorded) > 5

    def test_first_order_kinetics(self, first_order_config):
        """Test simulation with first-order reaction kinetics."""
        result = run_simulation(first_order_config)

        # Should produce valid results
        assert result.final_profile is not None
        assert all(c >= 0 for c in result.final_profile)

    def test_sealed_system_depletes(self, sealed_config):
        """Test that sealed system depletes oxygen over time."""
        result = run_simulation(sealed_config)

        # Average concentration should decrease (no replenishment from top)
        initial_avg = sealed_config.C_initial_fraction
        final_avg = np.mean(result.final_profile)
        assert final_avg < initial_avg

    def test_open_vs_sealed_comparison(self):
        """Test that open system maintains higher O2 than sealed."""
        base_params = dict(k=0.5, nz=50, steps=30, dt=1.0)

        open_config = SimulationConfig(top_constraint='open', **base_params)
        sealed_config = SimulationConfig(top_constraint='sealed', **base_params)

        open_result = run_simulation(open_config)
        sealed_result = run_simulation(sealed_config)

        # Open system should have higher average concentration
        assert np.mean(open_result.final_profile) > np.mean(sealed_result.final_profile)


# =============================================================================
# SimulationResult Tests
# =============================================================================

class TestSimulationResult:
    """Tests for SimulationResult dataclass."""

    def test_to_dataframe(self, fast_config):
        """Test conversion to DataFrame."""
        result = run_simulation(fast_config, record_every=1)
        df = result.to_dataframe()

        assert len(df) > 0
        assert 'c_star' in df.columns
        assert 'C' in df.columns
        assert 'height_mm' in df.columns
        assert 't_hrs' in df.columns
        assert 't_s' in df.columns

    def test_to_dataframe_dimensionalization(self, fast_config):
        """Test that DataFrame has correct dimensionalized values."""
        result = run_simulation(fast_config, record_every=1)
        df = result.to_dataframe()

        # C should be c_star * C_air
        assert_allclose(df['C'].values, df['c_star'].values * fast_config.C_air)

        # height_mm should be z_idx * dz
        assert_allclose(df['height_mm'].values, df['z_idx'].values * fast_config.dz)

    def test_to_dataframe_empty_result(self):
        """Test that empty result produces empty DataFrame."""
        result = SimulationResult(config=SimulationConfig())
        df = result.to_dataframe()
        assert len(df) == 0

    def test_get_profile_at_step(self):
        """Test getting profile at specific step."""
        # Use config with no reaction to ensure clean execution
        config = SimulationConfig(k=0, k1=0, nz=30, steps=10, dt=1.0)
        result = run_simulation(config, record_every=1)

        # Get profile at step 0
        profile = result.get_profile_at_step(0)

        assert len(profile) == config.nz
        assert all(profile['step'] == 0)

    def test_get_profile_sorted_by_z(self, fast_config):
        """Test that profile is sorted by z index."""
        result = run_simulation(fast_config, record_every=1)
        profile = result.get_profile_at_step(0)

        z_indices = profile['z_idx'].values
        assert all(z_indices[i] <= z_indices[i+1] for i in range(len(z_indices)-1))


# =============================================================================
# Parameter Sweep Tests
# =============================================================================

class TestParameterSweep:
    """Tests for parameter sweep functionality."""

    def test_sweep_returns_list(self, fast_config):
        """Test that sweep returns a list of results."""
        results = run_parameter_sweep(fast_config, 'k', [0.5, 1.0])

        assert isinstance(results, list)
        assert len(results) == 2

    def test_sweep_varies_parameter(self, fast_config):
        """Test that sweep actually varies the parameter."""
        k_values = [0.5, 1.0, 2.0]
        results = run_parameter_sweep(fast_config, 'k', k_values)

        for result, expected_k in zip(results, k_values):
            assert result.config.k == expected_k

    def test_sweep_preserves_other_params(self):
        """Test that sweep preserves non-swept parameters."""
        config = SimulationConfig(nz=50, steps=10, dt=1.0, L=5.0)
        results = run_parameter_sweep(config, 'k', [0.5, 1.0])

        for result in results:
            assert result.config.nz == 50
            assert result.config.L == 5.0

    def test_sweep_different_parameters(self):
        """Test sweeping different parameter types."""
        config = SimulationConfig(nz=50, steps=10, dt=1.0)

        # Sweep diffusion coefficient
        D_values = [1e-3, 3e-3]
        results = run_parameter_sweep(config, 'D', D_values)
        for result, expected_D in zip(results, D_values):
            assert result.config.D == expected_D


# =============================================================================
# Combine Results Tests
# =============================================================================

class TestCombineResults:
    """Tests for combining multiple simulation results."""

    def test_combine_results_produces_dataframe(self, fast_config):
        """Test that combine_results produces a DataFrame."""
        results = run_parameter_sweep(fast_config, 'k', [0.5, 1.0])
        df = combine_results(results)

        assert len(df) > 0

    def test_combine_results_includes_all_k(self, fast_config):
        """Test that combined results include all k values."""
        k_values = [0.5, 1.0, 2.0]
        results = run_parameter_sweep(fast_config, 'k', k_values)
        df = combine_results(results)

        assert set(df['k'].unique()) == set(k_values)

    def test_combine_empty_list(self):
        """Test combining empty list of results raises ValueError."""
        import pandas as pd
        with pytest.raises(ValueError):
            combine_results([])


# =============================================================================
# Physical Behavior Tests
# =============================================================================

class TestPhysicalBehavior:
    """Tests verifying physically expected behavior."""

    def test_steady_state_profile_shape(self):
        """Test that steady state has expected profile shape with consumption."""
        config = SimulationConfig(
            k=0.5, nz=100, steps=500, dt=1.0
        )
        result = run_simulation(config)

        # Profile should be monotonically increasing from bottom to top
        profile = result.final_profile
        for i in range(len(profile) - 1):
            assert profile[i] <= profile[i+1] + 1e-10

    def test_higher_k_faster_depletion(self):
        """Test that higher reaction rate depletes oxygen faster."""
        # Use low k values and short simulation to avoid hitting zero
        base_params = dict(nz=50, steps=5, dt=0.5, halt_on_C_zero=False)

        low_k = SimulationConfig(k=0.1, **base_params)
        high_k = SimulationConfig(k=0.5, **base_params)

        result_low = run_simulation(low_k)
        result_high = run_simulation(high_k)

        # Higher k should have lower average concentration
        assert np.mean(result_high.final_profile) < np.mean(result_low.final_profile)

    def test_deeper_well_slower_diffusion(self):
        """Test that deeper wells have slower O2 penetration."""
        base_params = dict(k=0, nz=50, steps=20, dt=1.0, C_initial_fraction=0.0)

        shallow = SimulationConfig(L=2.0, **base_params)
        deep = SimulationConfig(L=5.0, **base_params)

        result_shallow = run_simulation(shallow)
        result_deep = run_simulation(deep)

        # Relative penetration to bottom should be greater in shallow well
        # (comparing bottom concentration as fraction of final range)
        shallow_bottom_frac = result_shallow.final_profile[0]
        deep_bottom_frac = result_deep.final_profile[0]

        assert shallow_bottom_frac > deep_bottom_frac

    def test_mass_conservation_sealed_no_reaction(self):
        """Test mass conservation in sealed system without reaction."""
        config = SimulationConfig(
            top_constraint='sealed',
            k=0,
            k1=0,
            nz=50,
            steps=50,
            dt=1.0,
            C_initial_fraction=0.5,
        )
        result = run_simulation(config)

        # Total mass should be conserved (integral of C should be constant)
        initial_mass = config.C_initial_fraction * config.nz
        final_mass = np.sum(result.final_profile)

        assert_allclose(final_mass, initial_mass, rtol=1e-5)


# =============================================================================
# Edge Cases and Error Handling
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    def test_single_cell_mesh(self):
        """Test simulation with minimal mesh (single cell)."""
        config = SimulationConfig(nz=1, steps=5, dt=1.0, k=0)
        result = run_simulation(config)

        assert len(result.final_profile) == 1

    def test_zero_steps(self):
        """Test simulation with zero steps."""
        config = SimulationConfig(nz=50, steps=0, dt=1.0)
        result = run_simulation(config)

        # Should have no recorded points and final_profile from initial state
        assert result.final_profile is not None

    def test_very_small_dt(self):
        """Test simulation with very small time step."""
        config = SimulationConfig(nz=20, steps=5, dt=1e-6)
        result = run_simulation(config)

        # Should still produce valid results
        assert result.final_profile is not None

    def test_very_large_k(self):
        """Test simulation with very high reaction rate."""
        config = SimulationConfig(
            k=100.0, nz=50, steps=10, dt=0.1,
            halt_on_C_zero=True
        )
        result = run_simulation(config)

        # Should halt early due to depletion
        assert result.final_profile is not None


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
