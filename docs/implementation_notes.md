# Implementation Notes

The first implementation is intentionally small and local-machine friendly.

## Expected Raw Genotype Format

Consumer genotype text files are expected to contain at least:

```text
rsid chromosome position genotype
```

Comment lines starting with `#` are ignored. Common 23andMe-like tab-separated files are supported.

## Phenotype Format

The cleaned phenotype file should contain:

```text
sample_id,celiac_status
sample_001,case
sample_002,control
```

`celiac_status` must be `case` or `control`.

## Research-Only Warning

The project produces a research risk estimate only. It is not a diagnostic system.

