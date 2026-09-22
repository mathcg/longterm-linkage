# Reproduction materials

`original/` preserves the files distributed with the JASA article. The package
API in `R/api.R` provides a safer and more consistent interface to the same
estimators and data generators.

The source bundle was downloaded from the article software page maintained by
Ying-Qi Zhao:

<https://research.fredhutch.org/content/dam/research/zhao/publications/_codeLinkageJASASubmission.zip>

The SHA-256 checksum of the downloaded ZIP is
`e0febed8c48eb4fb649fa9f7ac31866df50a2dbc4f8b1aa9b8017cef995ebd9e`.

The original bundle contains:

- `linkage_functions.R`: simulation designs and estimators;
- `linkage_example.R`: repeated simulation driver;
- `simulation_result_collection.R`: simulation summaries;
- `run_all.R`: full simulation entry point;
- `SWOGData_analysis_functions.R`: PCPT analysis functions;
- `run_SWOGData_analysis.R`: synthetic PCPT analysis entry point; and
- two Word documents describing the simulation and PCPT pipelines.

The synthetic PCPT dataset is also installed with the package. Locate it with:

```r
system.file("extdata", "PCPT_syntheticdata.txt", package = "longtermLinkage")
```

Do not place original participant-level PCPT or Medicare data in this
repository.
