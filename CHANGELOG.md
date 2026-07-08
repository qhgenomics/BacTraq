# Changelog

All notable changes to BacTraq are documented in this file.

## [1.1.0-rc.3] - 2026-07-08

### Changed

- `--summary` now writes `<output>_summary.xlsx` instead of a `.csv`, with two sheets: `Events` (the existing per-lineage Status/New Name table) and `Rename Trace` (one row per sample, Old/New cluster name per threshold, coarsest first; singletons marked `Unclustered`, samples absent from `--history` marked `New Sample` in the Old column).
- Added `openpyxl` as a dependency (required for Excel output).

## [1.1.0-rc.2] - 2026-07-02

### Changed

- Refactored the monolithic `bactraq.py` into separate modules: `cli.py`, `clustering.py`, `cluster_trace.py`, `reporting.py`, and `helper.py` for clearer separation of concerns.
- Updated CLI help messages for clarity.
- Updated README and docs, including a new cluster-trace diagram (`docs/cluster-trace-diagram.svg`) replacing older PNG screenshots.
- Updated example run outputs and tests to match the refactored code.

### Fixed

- Fixed an error where a cascading parent cluster name was not updated correctly when a name changed.

## [1.1.0-rc.1] - 2026-07-01

### Added

- New logic for handling singleton clusters and logging cluster-naming events (merge/split/new/consistent decisions are now logged).
- Added pytest coverage for core clustering/naming functions.

### Changed

- Removed the pre-built wheel from version control (build artifacts no longer committed).

## [1.0.2] - 2026-01-20

### Added

- `--nRef` flag to support SNP-dist matrices without a `Reference` row/column.

### Changed

- Migrated packaging from `setup.py` to `pyproject.toml`.

## [1.0.1] and earlier

### Fixed

- Various fixes to history file handling and naming consistency across runs.

### Changed

- Documentation updates covering installation, in-house naming structure, and clustering naming rules.

[1.1.0-rc.2]: https://github.com/qhgenomics/BacTraq/compare/6b6d435...7628c08
[1.1.0-rc.1]: https://github.com/qhgenomics/BacTraq/compare/94edd38...6b6d435
[1.0.2]: https://github.com/qhgenomics/BacTraq/compare/96a6748...94edd38
