# HDStIM-paper
Code to reproduce analysis for the HDStIM paper.

## Bone marrow dataset
### HDStIM run and marker ranking
```
scripts/bone_marrow/HDStIM_run.R
scripts/bone_marrow/HDStIM_marker_ranking.R
```
### Figures from bone marrow dataset
#### Figures 1 B & E. Supplementary figures 1 A & B
```
scripts/bone_marrow/HDStIM_vs_Bendall_paper.R
```
#### Figure 1 C & D
```
scripts/bone_marrow/HDStIM_diagnostic_plots.R
```
#### Figure 1 F
```
scripts/bone_marrow/marker_ranking_z_score_heatmap.R
```
#### Supplementary figure 1 C
```
scripts/bone_marrow/stimulated_vs_unstimulated_fc_IFNa.R
```
#### Supplementary figure 1 D
```
scripts/bone_marrow/stimulated_well_heatmap_IFNa.R
```

## Pediatrics dataset
## HDStIM run and diagnostic plots (not in the paper)
```
scripts/pediatrics/HDStIM_run.R 
scripts/pediatrics/HDStIM_diagnostic_plots.R
scripts/pediatrics/corrections_to_HDStIM_output.R
```

### Figures from pediatrics dataset
#### Figure 2 A
```
scripts/pediatrics/pediatrics_age_distribution.R
```

#### Figure 2 B
```
scripts/pediatrics/resp_noresp_marker_ranking_heatmap.R
```

#### Figure 2 C
```
scripts/pediatrics/resp_delta_lm_frac_bubble.R
```

#### Figure 2 D
```
scripts/pediatrics/resp_unstim_spline_heatmap.R
```

#### Figure 2 E
```
scripts/pediatrics/resp_unstim_mag_ranking.R
```

#### Supplementary Figure 2 A
```
scripts/pediatrics/lm_frac_case_and_control.R
scripts/pediatrics/lm_frac_case_and_control_heatmaps.R
```

#### Supplementary Figure 2 B
```
scripts/pediatrics/lm_frac_case_v_control_heatmaps.R
```

#### Supplementary Figure 2 C
```
scripts/pediatrics/pre_HDStIM_UMAP.R
```
