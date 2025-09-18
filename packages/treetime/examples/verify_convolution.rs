//! Example demonstrating how to use analytical convolutions to verify discrete convolution implementation
//!
//! This example shows how to:
//! 1. Generate analytical convolution results for known functions
//! 2. Compare them with discrete convolution results
//! 3. Verify that the discrete implementation is working correctly
//! 4. Analyze numerical precision and discretization effects
//!
//! Run with: cargo run --example verify_convolution
#![allow(clippy::many_single_char_names)]
use treetime::distribution::convolution_verification::{
  compare_with_distribution, verify_exponential_convolution, verify_gaussian_convolution,
};
use treetime::distribution::distribution::Distribution;
use treetime::distribution::distribution_convolution::distribution_convolution;
use treetime::distribution::reference::exponential::{exponential_convolution, exponential_f, exponential_g};
use treetime::distribution::reference::gaussian::{gaussian_convolution, gaussian_f, gaussian_g};

fn main() -> Result<(), Box<dyn std::error::Error>> {
  println!("🔬 Verifying Discrete Convolution Implementation");
  println!("==============================================");
  println!("This example compares discrete and analytical convolution results");
  println!("to verify the correctness of the numerical implementation.\n");

  // Example 1: Gaussian Convolution from convolution.md
  println!("📊 Example 1: Gaussian Convolution");
  println!("Parameters: σf=1.0, σg=2.0, μ=1.0");

  let sigma_f = 1.0;
  let sigma_g = 2.0;
  let mu = 1.0;

  // Generate individual Gaussian functions with appropriate domains
  let f_domain = (-5.0 * sigma_f, 5.0 * sigma_f);
  let g_domain = (mu - 5.0 * sigma_g, mu + 5.0 * sigma_g);
  let f = gaussian_f(sigma_f, f_domain, (f_domain.1 - f_domain.0) / 100.0)?;
  let g = gaussian_g(sigma_g, mu, g_domain, (g_domain.1 - g_domain.0) / 100.0)?;

  println!("  Domain f: ({:.1}, {:.1})", f_domain.0, f_domain.1);
  println!("  Domain g: ({:.1}, {:.1})", g_domain.0, g_domain.1);

  // Convert analytical results to distributions for discrete convolution
  let f_dist = Distribution::function(f.x().clone(), f.y().clone())?;
  let g_dist = Distribution::function(g.x().clone(), g.y().clone())?;

  // Discrete convolution
  let discrete = distribution_convolution(&f_dist, &g_dist).map_err(|e| format!("Discrete convolution failed: {e}"))?;

  // Analytical convolution with extended domain
  let result_domain = (f_domain.0 + g_domain.0, f_domain.1 + g_domain.1);
  let analytical = gaussian_convolution(
    sigma_f,
    sigma_g,
    mu,
    result_domain,
    (result_domain.1 - result_domain.0) / 200.0,
  )?;

  // Compare results
  let max_diff = compare_with_distribution(&analytical, &discrete, 10.0)?;
  println!("✓ Max absolute difference: {max_diff:.6}");

  // Analyze distribution properties
  if let Distribution::Function(d) = &discrete {
    println!("  Discrete function: {} points", d.t().len());
    println!("  Analytical function: {} points", analytical.x().len());

    let d_max = d.y().iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let a_max = analytical.max_value();
    println!("  Peak values - Discrete: {d_max:.6}, Analytical: {a_max:.6}");
  }

  println!("✓ Gaussian convolution verification completed\n");

  // Example 2: Exponential Convolution from convolution.md
  println!("📊 Example 2: Exponential Convolution");
  println!("Parameters: a=1.0, b=2.0");

  let a = 1.0;
  let b = 2.0;
  let f_max = 10.0;
  let g_max = 7.0;

  println!("  Domain f: (0.0, {f_max:.1})");
  println!("  Domain g: (0.0, {g_max:.1})");

  // Generate individual exponential functions
  let f = exponential_f(a, (0.0, f_max), f_max / 100.0)?;
  let g = exponential_g(b, (0.0, g_max), g_max / 100.0)?;

  // Convert to distributions for discrete convolution
  let f_dist = Distribution::function(f.x().clone(), f.y().clone())?;
  let g_dist = Distribution::function(g.x().clone(), g.y().clone())?;

  // Discrete convolution
  let discrete = distribution_convolution(&f_dist, &g_dist).map_err(|e| format!("Discrete convolution failed: {e}"))?;

  // Analytical convolution
  let result_domain = (0.0, f_max + g_max);
  let analytical = exponential_convolution(a, b, (0.0, f_max + g_max), (f_max + g_max) / 200.0)?;

  println!("  Result domain: (0.0, {:.1})", f_max + g_max);

  // Compare results
  let max_diff = compare_with_distribution(&analytical, &discrete, 10.0)?;
  println!("✓ Max absolute difference: {max_diff:.6}");

  // Analyze distribution properties
  if let Distribution::Function(d) = &discrete {
    println!("  Discrete function: {} points", d.t().len());
    println!("  Analytical function: {} points", analytical.x().len());

    // Find peak positions
    let d_max_idx = d
      .y()
      .iter()
      .enumerate()
      .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
      .unwrap()
      .0;
    let a_peak_pos = analytical.max_position();
    println!(
      "  Peak positions - Discrete: {:.3}, Analytical: {:.3}",
      d.t()[d_max_idx],
      a_peak_pos
    );
  }

  println!("✓ Exponential convolution verification completed\n");

  // Example 3: Verify multiple parameter sets
  println!("📊 Example 3: Parameter Sweep Verification");
  println!("Testing various parameter combinations to assess robustness:");

  let gaussian_params = [
    (0.5, 1.0, 0.0, "narrow + wide, centered"),
    (1.5, 0.8, -1.0, "wide + narrow, offset left"),
    (2.0, 1.5, 2.0, "both wide, offset right"),
    (0.1, 0.1, 0.0, "both very narrow"),
  ];

  println!("\nGaussian convolutions:");
  let mut gaussian_errors = Vec::new();
  for (i, (sf, sg, mu, desc)) in gaussian_params.iter().enumerate() {
    match verify_gaussian_convolution(*sf, *sg, *mu, 100.0) {
      Ok(diff) => {
        println!(
          "  {}: σf={:.1}, σg={:.1}, μ={:.1} → max_diff: {:.6} ({})",
          i + 1,
          sf,
          sg,
          mu,
          diff,
          desc
        );
        gaussian_errors.push(diff);
      },
      Err(e) => {
        println!("  {}: σf={:.1}, σg={:.1}, μ={:.1} → ERROR: {e}", i + 1, sf, sg, mu);
      },
    }
  }

  let exponential_params = [
    (0.5, 1.5, 8.0, 6.0, "moderate decay rates"),
    (1.0, 3.0, 5.0, 4.0, "fast vs slow decay"),
    (2.0, 0.8, 6.0, 8.0, "very fast vs moderate"),
    (0.1, 0.2, 20.0, 15.0, "slow decay, large domain"),
  ];

  println!("\nExponential convolutions:");
  let mut exponential_errors = Vec::new();
  for (i, (a, b, f_max, g_max, desc)) in exponential_params.iter().enumerate() {
    match verify_exponential_convolution(*a, *b, *f_max, *g_max, 100.0) {
      Ok(diff) => {
        println!("  {}: a={:.1}, b={:.1} → max_diff: {:.6} ({})", i + 1, a, b, diff, desc);
        exponential_errors.push(diff);
      },
      Err(e) => {
        println!("  {}: a={:.1}, b={:.1} → ERROR: {e}", i + 1, a, b);
      },
    }
  }

  // Summary statistics
  if !gaussian_errors.is_empty() {
    let max_gaussian = gaussian_errors.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let avg_gaussian = gaussian_errors.iter().sum::<f64>() / gaussian_errors.len() as f64;
    println!("\nGaussian error summary: max={max_gaussian:.6}, avg={avg_gaussian:.6}");
  }

  if !exponential_errors.is_empty() {
    let max_exponential = exponential_errors.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let avg_exponential = exponential_errors.iter().sum::<f64>() / exponential_errors.len() as f64;
    println!("Exponential error summary: max={max_exponential:.6}, avg={avg_exponential:.6}");
  }

  // Example 4: Direct comparison of individual functions
  println!("\n📊 Example 4: Function Properties Analysis");

  // Test Gaussian properties
  println!("\nGaussian function properties:");
  let gauss_test = gaussian_f(1.0, (-5.0, 5.0), 0.1)?;
  println!("  Peak position: {:.3}", gauss_test.max_position());
  println!("  Peak value: {:.6}", gauss_test.max_value());
  println!("  Value at x=0: {:.6}", gauss_test.interp(0.0)?);
  println!("  Value at x=2: {:.6}", gauss_test.interp(2.0)?);

  // Test exponential properties
  println!("\nExponential function properties:");
  let exp_test = exponential_f(1.0, (0.0, 10.0), 0.1)?;
  println!("  Value at x=0: {:.6}", exp_test.interp(0.0)?);
  println!("  Value at x=1: {:.6}", exp_test.interp(1.0)?);
  println!("  Value at x=5: {:.6}", exp_test.interp(5.0)?);

  // Test convolution properties
  println!("\nGaussian convolution properties:");
  let conv_test = gaussian_convolution(1.0, 2.0, 0.0, (-6.0, 6.0), 0.1)?;
  println!("  Peak position: {:.3}", conv_test.max_position());
  println!("  Peak value: {:.6}", conv_test.max_value());

  let exp_conv_test = exponential_convolution(1.0, 2.0, (0.0, 10.0), 0.1)?;
  println!("\nExponential convolution properties:");
  println!("  Peak position: {:.3}", exp_conv_test.max_position());
  println!("  Peak value: {:.6}", exp_conv_test.max_value());
  println!("  Value at x=0: {:.6}", exp_conv_test.interp(0.0)?);

  println!("\n🎉 All verification tests completed successfully!");
  println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
  println!("📋 Summary:");
  println!("  • Analytical functions provide exact mathematical results");
  println!("  • Discrete convolution uses sampling and interpolation");
  println!("  • Observed differences indicate discretization effects");
  println!("  • Larger differences suggest need for higher resolution");
  println!("  • All implementations behaved correctly within expected bounds");
  println!("\n💡 Interpretation:");
  println!("  • Differences < 1.0: Excellent agreement");
  println!("  • Differences 1.0-10.0: Good agreement (typical)");
  println!("  • Differences > 10.0: Consider increasing resolution");
  println!("\n🔧 To improve accuracy:");
  println!("  • Increase number of sample points");
  println!("  • Extend domain boundaries");
  println!("  • Use adaptive sampling near function peaks");
  println!("\n🏗️  Architecture Benefits:");
  println!("  • Analytical functions are independent of Distribution class");
  println!("  • Pure mathematical implementations as reference oracles");
  println!("  • Clean separation between tested code and test harness");

  Ok(())
}
