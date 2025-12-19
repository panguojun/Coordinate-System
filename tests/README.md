# Tests

This directory contains verification tests for the Coordinate System Algebra framework.

## Running Tests

To run the geometry verification tests:

```bash
python tests/geometry_verification.py
```

## Test Coverage

The `geometry_verification.py` module provides comprehensive validation of:

1. **Classical Surface Curvature**
   - Sphere (K = 1/R²)
   - Cylinder (K = 0)
   - Hyperboloid (saddle, K < 0)
   - Torus (variable K)

2. **Accuracy Validation**
   - Machine precision: relative error < 10⁻⁹ to 10⁻¹⁷
   - Rotational invariance
   - Scale covariance

3. **Performance Testing**
   - Computational complexity (O(1) per point)
   - Comparison with traditional methods

## Expected Output

All tests should pass with:
- Pass rate: 100%
- Relative errors < 10⁻⁶ for all surfaces
- Standard deviation < 10⁻¹⁰ for symmetry tests

## Requirements

```bash
pip install numpy
```
