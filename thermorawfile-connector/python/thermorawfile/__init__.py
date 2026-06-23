"""Pure-Rust reader and writer for Thermo Finnigan ``.raw`` mass-spectrometry files.

No .NET runtime, no ``RawFileReader`` DLL — a clean-room implementation. The heavy
lifting lives in the compiled ``_core`` extension; this module is the Pythonic wrapper
(typed classes, properties, context manager, numpy peaks).

    >>> import thermorawfile as trf
    >>> rf = trf.RawFile("run.raw")
    >>> rf.instrument_model, rf.acquired
    ('Orbitrap Astral', '2024-11-03 13:09:44')
    >>> mz, intensity = rf.peaks(2)          # numpy arrays
    >>> rf.scan_filter(2)
    'FTMS + p NSI Full ms [350.00-1500.00]'
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterator, NamedTuple, Optional, Sequence, Union

import numpy as np
from numpy.typing import NDArray

from ._core import RawFile as _RawFile
from ._core import __version__ as __version__

__all__ = ["RawFile", "ScanEvent", "Peaks", "__version__"]

_Floats = Union[Sequence[float], NDArray]


class Peaks(NamedTuple):
    """A centroid peak list. Unpacks as ``(mz, intensity)``; also ``.mz`` / ``.intensity``."""

    mz: NDArray[np.float64]
    intensity: NDArray[np.float32]


@dataclass(frozen=True)
class ScanEvent:
    """Scan-event metadata for one scan."""

    ms_order: int
    analyzer: int
    isolation_center: float
    isolation_width: float
    collision_energy: float


class RawFile:
    """A Thermo Finnigan ``.raw`` file (revision >= 64, Orbitrap-era).

    The whole file is held in memory; reads are per scan. Writing is *functional* —
    replace or overlay a scan's centroids, then :meth:`save`. There are deliberately no
    setters for forensic provenance fields (serial number, audit timestamps): a writable
    open format cannot enforce authenticity, so that is a job for cryptographic provenance,
    not this writer.
    """

    __slots__ = ("_inner",)

    def __init__(self, path: str) -> None:
        self._inner = _RawFile(path)

    # -- dunder --
    def __len__(self) -> int:
        return self._inner.n_scans

    def __repr__(self) -> str:
        return repr(self._inner)

    def __iter__(self) -> Iterator[int]:
        """Iterate scan numbers ``first_scan`` .. ``last_scan`` (lightweight — no spectra)."""
        return iter(range(self._inner.first_scan, self._inner.last_scan + 1))

    def __enter__(self) -> "RawFile":
        return self

    def __exit__(self, *exc: object) -> None:
        return None

    # -- properties --
    @property
    def version(self) -> int:
        """File-format revision (e.g. 66)."""
        return self._inner.version

    @property
    def n_scans(self) -> int:
        return self._inner.n_scans

    @property
    def first_scan(self) -> int:
        return self._inner.first_scan

    @property
    def last_scan(self) -> int:
        return self._inner.last_scan

    @property
    def instrument_model(self) -> Optional[str]:
        """Detected instrument model, e.g. ``"Orbitrap Astral"`` (or ``None``)."""
        return self._inner.instrument_model

    @property
    def acquired(self) -> Optional[str]:
        """Acquisition date ``"YYYY-MM-DD HH:MM:SS"`` (or ``None``)."""
        return self._inner.acquired

    @property
    def path(self) -> str:
        return self._inner.path

    # -- reads --
    def checksum_valid(self) -> bool:
        """Whether the file's (keyless) Adler-32 integrity checksum matches its content."""
        return self._inner.checksum_valid()

    def peaks(self, scan: int) -> Peaks:
        """Centroid peaks for ``scan`` (1-based) as a :class:`Peaks` of numpy arrays."""
        mz, intensity = self._inner.peaks(scan)
        return Peaks(mz, intensity)

    def scan_filter(self, scan: int) -> Optional[str]:
        """The Thermo filter line, e.g. ``"FTMS + p NSI Full ms [350.00-1500.00]"``."""
        return self._inner.scan_filter(scan)

    def scan_event(self, scan: int) -> Optional[ScanEvent]:
        """Scan-event metadata as a :class:`ScanEvent` (or ``None``)."""
        d = self._inner.scan_event(scan)
        return ScanEvent(**d) if d is not None else None

    def scan_params(self, scan: int) -> Optional[Dict[str, object]]:
        """Per-scan trailer parameters as ``{label: value}`` (or ``None``).

        Labels are the raw Thermo strings, e.g. ``"Ion Injection Time (ms):"``,
        ``"Charge State:"``, ``"FT Resolution:"``, ``"FAIMS CV:"``.
        """
        return self._inner.scan_params(scan)

    # -- writes (functional only) --
    def author_centroids(self, scan: int, mz: _Floats, intensity: _Floats) -> None:
        """Replace ``scan``'s centroid peaks. ``mz``/``intensity`` may be lists or arrays
        (equal length; must fit the scan's existing packet budget)."""
        self._inner.author_centroids(
            scan,
            np.ascontiguousarray(mz, dtype=np.float64),
            np.ascontiguousarray(intensity, dtype=np.float32),
        )

    def overlay_centroids(
        self, scan: int, mz: _Floats, intensity: _Floats, merge_tol_ppm: float = 20.0
    ) -> None:
        """Merge synthetic peaks onto ``scan``'s existing centroids (within ``merge_tol_ppm``)."""
        self._inner.overlay_centroids(
            scan,
            np.ascontiguousarray(mz, dtype=np.float64),
            np.ascontiguousarray(intensity, dtype=np.float32),
            merge_tol_ppm,
        )

    def save(self, path: str) -> None:
        """Recompute the integrity checksum and write the file to ``path``."""
        self._inner.save(path)
