#[macro_use]
extern crate criterion;

use criterion::Criterion;
use rand::distributions::{Distribution, Uniform};

fn criterion_benchmark(c: &mut Criterion) {
    // f32
    for &n in &[100, 1000, 10000] {
        let in_ = {
            let mut rng = rand::thread_rng();
            let between = Uniform::from(0.0..2.0 * std::f32::consts::PI);
            let mut buf = vec![0.0; n];
            for val in buf.iter_mut() {
                *val = between.sample(&mut rng);
            }
            buf
        };

        let mut out = vec![0.0_f32; n];
        c.bench_function(&format!("cos32_n{}", n), |b| {
            b.iter(|| {
                for i in 0..n {
                    out[i] = in_[i].cos();
                }
            })
        });
        c.bench_function(&format!("vcos32_n{}", n), |b| {
            b.iter(|| unsafe {
                intel_mkl_sys::vsCos(n as i32, in_.as_ptr(), out.as_mut_ptr());
            })
        });
    }

    // f64
    for &n in &[100, 1000, 10000] {
        let in_ = {
            let mut rng = rand::thread_rng();
            let between = Uniform::from(0.0..2.0 * std::f64::consts::PI);
            let mut buf = vec![0.0; n];
            for val in buf.iter_mut() {
                *val = between.sample(&mut rng);
            }
            buf
        };

        let mut out = vec![0.0_f64; n];
        c.bench_function(&format!("cos64_n{}", n), |b| {
            b.iter(|| {
                for i in 0..n {
                    out[i] = in_[i].cos();
                }
            })
        });

        c.bench_function(&format!("vcos64_n{}", n), |b| {
            b.iter(|| unsafe {
                intel_mkl_sys::vdCos(n as i32, in_.as_ptr(), out.as_mut_ptr());
            })
        });
    }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
