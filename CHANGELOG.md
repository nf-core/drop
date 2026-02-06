# nf-core/drop: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of nf-core/drop, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`
- Fixed BAM index file renaming issue in `countReads.R` that caused "valid 'index' file required" errors ([#XXX](https://github.com/nf-core/drop/pull/XXX))
- Switched to SerialParam for BiocParallel to improve compatibility with containerized environments
### `Dependencies`

### `Deprecated`
