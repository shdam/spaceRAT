# spaceRAT 0.99.0

* Streamline plot code in `projectSample`
* Repel labels on loading plot
* Improve new sample coloring when projected
* Adjusted to updated `spaceRATScaffolds`
* Added `add_umap` option for `buildScaffold` to include UMAP scaffolds
* Parameter renaming: `space` -> `scaffold` and `dim_reduction` -> `dimred`
* Added parameter that enables or disables differential expression analysis in `buildScaffold` (`subset_deg`).
* Fixed `as.data.frame` bug where rownames are not assigned when converting from tibbles
* "logged" -> "exprs"
* Preparing package for Bioconductor submission
* Wrote vignettes
* Added tests
* Moved scaffolds and example data to a separate data package called 
[spaceRATScaffolds](https://github.com/shdam/spaceRATScaffolds)
