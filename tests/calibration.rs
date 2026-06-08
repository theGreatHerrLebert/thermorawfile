//! Frequency↔m/z calibration round-trip (no data file needed).
//!
//! Uses the real coefficients decoded from a rev66 Orbitrap Astral MS1 scan
//! (PRIDE PXD074900): nparam=5, m/z = b/f² + c/f⁴.

use thermorawfile::Calibration;

#[test]
fn nparam5_inverse_round_trips() {
    let cal = Calibration {
        nparam: 5,
        a: 0.0,
        b: 1.6955368e8,
        c: 1.3376782e8,
    };
    for &mz in &[401.0, 500.0, 533.1506647, 650.0, 800.0, 899.0] {
        let f = cal.freq(mz).expect("reachable");
        let back = cal.mz(f);
        assert!(
            (back - mz).abs() < 1e-6,
            "m/z {mz} -> f {f} -> m/z {back} (Δ {})",
            (back - mz).abs()
        );
    }
}

#[test]
fn nparam4_inverse_round_trips() {
    let cal = Calibration {
        nparam: 4,
        a: -0.004,
        b: 3.34,
        c: 1.6955e8,
    };
    for &mz in &[300.0, 600.0, 1200.0] {
        let f = cal.freq(mz).expect("reachable");
        assert!((cal.mz(f) - mz).abs() < 1e-6);
    }
}
