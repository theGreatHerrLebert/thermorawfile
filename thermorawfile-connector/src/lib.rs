//! Python (PyO3) bindings for thermorawfile. Exposes the pure-Rust reader as
//! `thermorawfile.RawFile`, with numpy-backed peaks and dict-shaped metadata.
use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1};
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use thermorawfile::generic_record::GenericValue;
use thermorawfile::RawFile;

/// A parsed Thermo Finnigan `.raw` file.
#[pyclass(name = "RawFile")]
pub struct PyRawFile {
    inner: RawFile,
    path: String,
}

#[pymethods]
impl PyRawFile {
    #[new]
    fn new(py: Python<'_>, path: String) -> PyResult<Self> {
        let p = path.clone();
        let inner = py
            .allow_threads(|| RawFile::open(&p))
            .map_err(|e| PyIOError::new_err(e.to_string()))?;
        Ok(Self { inner, path })
    }

    #[getter]
    fn version(&self) -> u32 {
        self.inner.version
    }
    #[getter]
    fn n_scans(&self) -> usize {
        self.inner.scan_count()
    }
    #[getter]
    fn first_scan(&self) -> u32 {
        self.inner.first_scan
    }
    #[getter]
    fn last_scan(&self) -> u32 {
        self.inner.last_scan
    }
    #[getter]
    fn instrument_model(&self) -> Option<String> {
        self.inner.instrument_model.map(str::to_string)
    }
    #[getter]
    fn acquired(&self) -> Option<String> {
        self.inner.acquired.map(|d| d.to_iso())
    }
    #[getter]
    fn path(&self) -> &str {
        &self.path
    }

    fn __len__(&self) -> usize {
        self.inner.scan_count()
    }
    fn __repr__(&self) -> String {
        format!(
            "RawFile(rev {}, {} scans, {})",
            self.inner.version,
            self.inner.scan_count(),
            self.inner.instrument_model.unwrap_or("unknown model"),
        )
    }

    /// Whether the file's Adler-32 integrity checksum matches its content.
    fn checksum_valid(&self) -> bool {
        self.inner.checksum_valid()
    }

    /// Centroid peaks for `scan` (1-based) as `(mz: float64[], intensity: float32[])`.
    fn peaks(&self, py: Python<'_>, scan: u32) -> (Py<PyArray1<f64>>, Py<PyArray1<f32>>) {
        let peaks = self.inner.centroid_peaks(scan);
        let mz: Vec<f64> = peaks.iter().map(|p| p.mz).collect();
        let intensity: Vec<f32> = peaks.iter().map(|p| p.intensity).collect();
        (mz.into_pyarray(py).unbind(), intensity.into_pyarray(py).unbind())
    }

    /// The human-readable scan filter line, e.g. `"FTMS + p NSI Full ms [350.00-1200.00]"`.
    fn scan_filter(&self, scan: u32) -> Option<String> {
        self.inner.scan_filter(scan)
    }

    /// Scan-event metadata dict: ms_order, analyzer, isolation center/width, collision energy.
    fn scan_event<'py>(&self, py: Python<'py>, scan: u32) -> PyResult<Option<Bound<'py, PyDict>>> {
        let Some(e) = self.inner.scan_event(scan) else {
            return Ok(None);
        };
        let d = PyDict::new(py);
        d.set_item("ms_order", e.ms_order)?;
        d.set_item("analyzer", e.analyzer)?;
        d.set_item("isolation_center", e.isolation_center)?;
        d.set_item("isolation_width", e.isolation_width)?;
        d.set_item("collision_energy", e.collision_energy)?;
        Ok(Some(d))
    }

    /// Per-scan trailer parameters as a `{label: value}` dict (AGC, ion-injection time,
    /// charge, FT resolution, FAIMS, NCE, ...). `None` if the scan has no trailer record.
    fn scan_params<'py>(&self, py: Python<'py>, scan: u32) -> PyResult<Option<Bound<'py, PyDict>>> {
        let Some(p) = self.inner.scan_params(scan) else {
            return Ok(None);
        };
        let d = PyDict::new(py);
        for (label, val) in &p.record().values {
            match val {
                GenericValue::Int8(x) => d.set_item(label, *x)?,
                GenericValue::UInt8(x) => d.set_item(label, *x)?,
                GenericValue::Int16(x) => d.set_item(label, *x)?,
                GenericValue::UInt16(x) => d.set_item(label, *x)?,
                GenericValue::Int32(x) => d.set_item(label, *x)?,
                GenericValue::UInt32(x) => d.set_item(label, *x)?,
                GenericValue::Float32(x) => d.set_item(label, *x)?,
                GenericValue::Float64(x) => d.set_item(label, *x)?,
                GenericValue::Bool(x) => d.set_item(label, *x)?,
                GenericValue::String(s) => d.set_item(label, s)?,
                GenericValue::Gap => {}
            }
        }
        Ok(Some(d))
    }

    // --- writer: functional only (no forensic-marker setters — project policy) ---

    /// Replace scan `scan`'s centroid peaks with `(mz: float64[], intensity: float32[])`.
    /// Arrays must be equal length, with finite positive m/z and non-negative intensity.
    fn author_centroids(
        &mut self,
        scan: u32,
        mz: PyReadonlyArray1<f64>,
        intensity: PyReadonlyArray1<f32>,
    ) -> PyResult<()> {
        let peaks = peak_pairs(&mz, &intensity)?;
        self.inner.author_centroids(scan, &peaks).map_err(io_err)
    }

    /// Merge synthetic peaks onto scan `scan`'s existing centroids (within `merge_tol_ppm`).
    #[pyo3(signature = (scan, mz, intensity, merge_tol_ppm = 20.0))]
    fn overlay_centroids(
        &mut self,
        scan: u32,
        mz: PyReadonlyArray1<f64>,
        intensity: PyReadonlyArray1<f32>,
        merge_tol_ppm: f64,
    ) -> PyResult<()> {
        let peaks = peak_pairs(&mz, &intensity)?;
        self.inner
            .overlay_centroids(scan, &peaks, merge_tol_ppm)
            .map_err(io_err)
    }

    /// Recompute the integrity checksum and write the (possibly modified) file to `path`.
    fn save(&mut self, path: String) -> PyResult<()> {
        self.inner.save(&path).map_err(io_err)
    }
}

fn io_err(e: std::io::Error) -> PyErr {
    PyIOError::new_err(e.to_string())
}

fn peak_pairs(
    mz: &PyReadonlyArray1<f64>,
    intensity: &PyReadonlyArray1<f32>,
) -> PyResult<Vec<(f64, f32)>> {
    let mz = mz.as_slice()?;
    let inten = intensity.as_slice()?;
    if mz.len() != inten.len() {
        return Err(PyValueError::new_err(format!(
            "mz ({}) and intensity ({}) must have the same length",
            mz.len(),
            inten.len()
        )));
    }
    if mz
        .iter()
        .zip(inten)
        .any(|(&m, &i)| !m.is_finite() || m <= 0.0 || !i.is_finite() || i < 0.0)
    {
        return Err(PyValueError::new_err(
            "peaks require finite positive m/z and finite non-negative intensity",
        ));
    }
    Ok(mz.iter().zip(inten).map(|(&m, &i)| (m, i)).collect())
}

/// The compiled extension, imported by the `thermorawfile` Python wrapper as `._core`.
#[pymodule]
fn _core(_py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyRawFile>()?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
