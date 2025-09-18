# Changelog

All notable changes to **GOTTA** will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

> Legend of sections:
> - **Added** for new features
> - **Changed** for behavior changes
> - **Deprecated** for soon-to-be removed
> - **Removed** for now-removed features
> - **Fixed** for any bug fixes
> - **Security** for vulnerabilities
> - **Docs** for documentation-only changes
> - **Performance** for speed/efficiency improvements
> - **Build/CI** for packaging, wheels, runners, etc.
> - **Breaking** for migration notes when APIs change (MAJOR bumps)

## [Unreleased]
### Added
- *(placeholder)*

### Changed
- *(placeholder)*

### Deprecated
- *(placeholder)*

### Removed
- *(placeholder)*

### Fixed
- *(placeholder)*

### Performance
- *(placeholder)*

### Docs
- *(placeholder)*

### Build/CI
- *(placeholder)*

### Security
- *(placeholder)*

### Known Issues
- *(placeholder — link to GitHub issues)*

### Breaking / Migration Guide
- *(placeholder — describe steps users must take to upgrade)*


---

## [0.1.0] - 2025-09-18
**Status:** Alpha (APIs may change)

### Added
- Initial project scaffold for **GOTTA** (Geometric Optimal Transport Tableau & Alignment):
  - Core primitives for **visualization**, **integration**, and **clustering** of spatially resolved omics.
  - Python package skeleton (`python/gotta/`) with stub APIs:
    - `gotta.visualize.*`, `gotta.integrate.*`, `gotta.cluster.*`
  - R bindings scaffold (`/R`) for Seurat-friendly data gateways.
  - Example datasets (toy) and notebooks in `/examples/`.
  - Basic docs site skeleton in `/docs/` with quickstart tutorial.
- Continuous Integration (CI) for linting & tests on Linux/macOS (GitHub Actions).
- MIT License and project governance docs (CODE_OF_CONDUCT, CONTRIBUTING).

### Docs
- Added `README` with installation notes, roadmap, and contribution guidelines.
- Added `CODE_OF_CONDUCT.md`, `CONTRIBUTING.md`, and `SECURITY.md`.

### Build/CI
- Python packaging metadata (`pyproject.toml`), test skeleton with `pytest`.
- R package scaffolding with `DESCRIPTION`, `NAMESPACE` (placeholders).
- Optional C++ toolchain checks in CI (cache build artifacts).

### Known Issues
- APIs are provisional; plotting themes and Seurat/AnnData adapters may change.

---

## Versioning Policy

- **MAJOR**: Breaking API changes or file format changes (e.g., `v1.0.0`).
- **MINOR**: Backward-compatible features and improvements (e.g., `v1.1.0`).
- **PATCH**: Backward-compatible bug fixes and small improvements (e.g., `v1.1.1`).

When **Breaking** changes occur, a **Migration Guide** will be included under each release with concrete steps (code snippets) to upgrade.

---

## Compare Links

- [Unreleased]: https://github.com/YOUR-ORG/GOTTA/compare/v0.1.0...HEAD
- [0.1.0]: https://github.com/YOUR-ORG/GOTTA/releases/tag/v0.1.0
