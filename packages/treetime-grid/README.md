# treetime-grid

Data structures and interpolation routines for representing and evaluating functions on uniformly-spaced grids.

## Types

- `Grid<T>` - Uniform grid parameters (x_min, dx, n_points) with indexing and iteration
- `GridFn<T>` - Function values on a uniform grid with piecewise linear interpolation
- `GridFnF64` - Type alias for `GridFn<f64>`
- `InterpElem` - Trait bound for grid element types (implemented for `f64`)
- `GridIter<T>` - Iterator over grid x coordinates

## Grid

Uniform grid parameters and coordinate calculations. Requires at least 2 points and positive spacing.

### Construction

```rust
use treetime_grid::Grid;
use ndarray::array;

let grid = Grid::from_start_dx(0.0, 0.1, 101)?;              // start, spacing, count
let grid = Grid::from_range_n_points(0.0, 10.0, 101)?;        // start, end, count
let grid = Grid::from_range_dx(0.0, 10.0, 0.1)?;              // start, end, spacing
let grid = Grid::from_array(&array![0.0, 0.5, 1.0, 1.5])?;   // from uniform array
```

### Properties and indexing

```rust
let x_min = grid.x_min();
let x_max = grid.x_max();
let (lo, hi) = grid.x_range();
let dx = grid.dx();
let n = grid.n_points();

let x = grid.x_at(50);                    // coordinate at index
let idx = grid.find_interval_index(5.5);   // interval containing value
let arr = grid.to_array();                 // all coordinates as Array1
```

### Iteration

```rust
for x in grid.iter() { }
for x in &grid { }
```

## GridFn

Function represented as values on a uniform grid with piecewise linear interpolation.

### Construction

```rust
use treetime_grid::{Grid, GridFn};
use ndarray::array;

// From a closure over a range
let f = GridFn::from_n_points((0.0, 1.0), 101, |x| x * x)?;
let f = GridFn::from_grid((0.0, 1.0), 0.01, |x| x.sin())?;

// From a Grid and closure
let grid = Grid::from_range_n_points(0.0, 1.0, 101)?;
let f = GridFn::from_grid_fn(grid, |x| x * x)?;

// From arrays
let f = GridFn::from_range_values((0.0, 2.0), array![0.0, 1.0, 4.0])?;
let f = GridFn::from_start_dx_values(0.0, 1.0, array![0.0, 1.0, 4.0])?;
let f = GridFn::from_grid_array(grid, array![/* ... */])?;
let f = GridFn::from_arrays(&x_uniform, y)?;

// Constants
let f = GridFn::constant((0.0, 1.0), 101, 5.0)?;
let f = GridFn::zeros((0.0, 1.0), 101)?;
let f = GridFn::ones((0.0, 1.0), 101)?;

// From non-uniform data (resamples to uniform grid)
let x = array![0.0, 0.1, 0.5, 1.0];
let y = array![0.0, 0.2, 1.0, 2.0];
let f = GridFn::from_arrays_nonuniform(&x, &y)?;
```

### Interpolation

```rust
let y = f.interp(0.5)?;
let ys = f.interp_many(&array![0.25, 0.5, 0.75])?;
```

### Resampling

```rust
let g = f.resample(&new_grid)?;
let g = f.resample_range_n_points((0.0, 1.0), 201)?;
let g = f.resample_range_dx((0.0, 1.0), 0.005)?;
let g = f.resample_start_dx(0.0, 0.005, 201)?;
```

### Properties

```rust
let grid = f.grid();
let y = f.y();                  // &Array1<T>
let x = f.x();                  // Array1<T> (allocates)
let (lo, hi) = f.x_range();
let (lo, hi) = f.y_range();
let pairs = f.to_pairs();       // Vec<(T, T)>
```

### Transformations

```rust
let g = f.mapv(|y| y * 2.0);          // transform values (new GridFn)
f.mapv_inplace(|y| y * 2.0);          // transform values in place
let g = f.negate_arg();               // f(x) -> f(-x), reflects across y-axis
f.negate_arg_inplace();               // reflect in place
```

## Interpolation behavior

- **Interior**: Piecewise linear interpolation between grid points
- **Extrapolation**: Constant (returns boundary values outside grid range)
- **Non-uniform resampling** (`interp_nonuniform`): Linear extrapolation outside input range
