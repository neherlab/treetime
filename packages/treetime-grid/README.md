# treetime-grid

Grid utilities for numerical computations on uniform grids. Provide data structures and interpolation routines for representing and evaluating functions on uniformly-spaced grids.

## Features

- **Grid**: Uniform grid parameters (x_min, dx, n_points) with indexing and iteration
- **GridFn**: Function values on a uniform grid with piecewise linear interpolation
- **Non-uniform resampling**: Convert non-uniformly spaced data to uniform grids

## API

### Grid

Encapsulate uniform grid parameters and provide coordinate calculations.

```rust
use treetime_grid::Grid;

// Create from start, spacing, and point count
let grid = Grid::from_start_dx(0.0, 0.1, 101)?;

// Create from range and point count
let grid = Grid::from_range_n_points(0.0, 10.0, 101)?;

// Create from range and spacing
let grid = Grid::from_range_dx(0.0, 10.0, 0.1)?;

// Access grid properties
let x_min = grid.x_min();
let x_max = grid.x_max();
let dx = grid.dx();
let n = grid.n_points();

// Get coordinate at index
let x = grid.x_at(50);

// Find interval containing a value
let idx = grid.find_interval_index(5.5);

// Iterate over coordinates
for x in grid.iter() {
    // ...
}
```

### GridFn

Represent a function as values on a uniform grid with linear interpolation.

```rust
use treetime_grid::GridFn;
use ndarray::array;

// Create from a function
let f = GridFn::from_n_points((0.0, 1.0), 101, |x| x * x)?;

// Create from arrays
let f = GridFn::from_range_values((0.0, 2.0), array![0.0, 1.0, 4.0])?;

// Interpolate at a point (linear interpolation, constant extrapolation)
let y = f.interp(0.5)?;

// Interpolate at multiple points
let ys = f.interp_many(&array![0.25, 0.5, 0.75])?;

// Resample to a different grid
let f_resampled = f.resample_range_n_points((0.0, 1.0), 201)?;

// Transform values
let f_scaled = f.mapv(|y| y * 2.0);

// Reflect across y-axis: f(x) -> f(-x)
let f_reflected = f.negate_arg();
```

### Non-uniform data

Convert non-uniformly spaced (x, y) pairs to a uniform grid.

```rust
use treetime_grid::GridFn;
use ndarray::array;

let x = array![0.0, 0.1, 0.5, 1.0];  // non-uniform spacing
let y = array![0.0, 0.2, 1.0, 2.0];

// Automatically resamples to uniform grid preserving detail
let f = GridFn::from_arrays_nonuniform(&x, &y)?;
```

## Interpolation behavior

- **Interior points**: Piecewise linear interpolation between grid points
- **Extrapolation**: Constant extrapolation (returns boundary values outside grid range)
