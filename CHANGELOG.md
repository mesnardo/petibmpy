# Change Log

---

## Current development

---

### Added

* Documentation with Mkdocs (built on gh-pages using doctr during Travis CI job).
* Code coverage with coveralls.
* Unit-tests
* Module to parse a PETSc log view file.
* Function to interpolate a 2D field from one grid to another.
* Module to create probes (volume and point) and to load data from files.
* Mesh feature to split the configuration of a stretched segment into two sub-configurations: one for a stretched segment where the maximum grid spacing is enforced, followed by one for a uniform segment (with a grid spacing equal to the maximum spacing of the previous stretched segment).
* `CartesianGrid` method to plot the gridlines of a 2D mesh in a Matplotlib figure.
* `CartesianGrid` method to print some information about the grid spacings of a mesh (min, max, ratio).

### Changed

* Update package version of dependencies.

### Fixed

### Removed

---

## 0.1

---

### Added

* Add modules with PEP8-compliant style.
* Add configuration file for Travis CI.

### Changed

### Fixed

### Removed
