# Changelog

All notable changes to sandlersteam will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.12.3] - 2026-02-21

### Fixed

- Bug fix in type conversion from superheated/subcooled inputs.

## [0.12.2] - 2026-02-07

### Fixed

- Included missing delta properties from saturated steam tables.

## [0.12.1] - 2026-02-04

### Fixed

- Bug in determining which table to use for supercritical conditions.
- Full incorporation of pint units throughout.

## [0.12.0] - 2026-02-04

### Changed

- CLI refactored for consistency.
- `State` class improvements.

## [0.11.2] - 2026-02-03

### Changed

- Package metadata update.

## [0.11.1] - 2026-02-03

### Changed

- `State` class and test suite updates.

## [0.10.0] - 2026-01-17

### Changed

- Offloaded thermodynamic state representation to `sandlermisc`; `State` now subclasses `ThermodynamicState` from `sandlermisc`.

## [0.9.0] - 2026-01-16

### Changed

- Reimplemented to use pint for units throughout.

## [0.8.1] - 2026-01-13

### Fixed

- Documentation configuration fix.

## [0.8.0] - 2026-01-13

### Changed

- Full reimplementation.

### Added

- ReadTheDocs documentation setup.

## [0.7.0] - 2025-12-31

### Changed

- Package metadata update.

## [0.6.1] - 2025-12-31

### Fixed

- CLI and saturated table fixes.

## [0.6.0] - 2025-12-31

### Changed

- CLI improvements; state handling updates.

## [0.5.1] - 2025-12-24

### Fixed

- Saturated and state module fixes.

## [0.5.0] - 2025-12-24

### Changed

- Migrated to `src/` layout; saturated steam data module restructured.

## [0.4.2] - 2025-03-17

### Fixed

- `State` and superheated table (`suph`) fixes.

## [0.4.1] - 2025-03-17

### Fixed

- `State` class and test suite updates.

## [0.4.0] - 2025-03-17

### Changed

- `State` class introduced.

## [0.3.3] - 2025-03-14

### Fixed

- Request module fixes.

## [0.3.2] - 2025-03-11

### Fixed

- Request module fixes.

## [0.3.1] - 2025-03-11

### Fixed

- Minor packaging fixes.

## [0.3.0] - 2025-03-11

### Added

- Initial packaged release.

## [0.1.0] - 2024-02-27

### Added

- Initial release of sandlersteam.
- Python API for programmatic use.