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
scripts/pediatric/HDStIM_run.R 
scripts/pediatric/HDStIM_diagnostic_plots.R
scripts/pediatric/corrections_to_HDStIM_output.R
```

### Figures from pediatric dataset
#### Figure 2 A
```
scripts/pediatric/pediatric_age_distribution.R
```

#### Figure 2 B
```
scripts/pediatric/resp_noresp_marker_ranking_heatmap.R
```

#### Figure 2 C
```
scripts/pediatric/resp_unstim_frac_lm_rank.R
scripts/pediatric/delta_frac_lm_rank.R
scripts/pediatric/resp_unstim_delta_lm_frac_bubble.R
```

#### Figure 2 D
```
scripts/pediatric/resp_unstim_spline_heatmap.R
```

#### Figure 2 E
```
scripts/pediatric/resp_unstim_mag_ranking.R
```

#### Supplementary Figure 2 A
```
scripts/pediatric/lm_frac_case_and_control.R
scripts/pediatric/lm_frac_case_and_control_heatmaps.R
```

#### Supplementary Figure 2 B
```
scripts/pediatric/lm_frac_case_v_control.R
scripts/pediatric/lm_frac_case_v_control_heatmaps.R
```

#### Supplementary Figure 2 C
```
scripts/pediatric/pre_HDStIM_UMAP.R
```
