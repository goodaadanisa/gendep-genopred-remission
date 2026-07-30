# Dissertation tables

Five main and fifteen supplementary tables from the final dissertation. Displayed values are rounded for readability; exact unrounded values remain in `table_data.yml`.

## Main tables

### Table 1. Baseline characteristics of the analysis cohort, overall and by treatment.

| Characteristic | Overall (n=430) | Treatment group 1 (n=210) | Treatment group 2 (n=220) |
| --- | --- | --- | --- |
| Age, years | 42.23 (11.43) | 42.77 (11.46) | 41.71 (11.41) |
| Age at onset, years | 31.38 (10.79) | 30.61 (10.58) | 32.11 (10.95) |
| Body mass index | 25.65 (4.81) | 25.81 (5.27) | 25.49 (4.32) |
| MADRS baseline score | 29.37 (6.56) | 29.83 (6.65) | 28.92 (6.45) |
| HRSD-17 baseline score | 24.30 (5.99) | 24.51 (6.26) | 24.11 (5.73) |
| BDI baseline score | 28.73 (9.29) | 28.90 (9.61) | 28.58 (8.99) |
| Sex code 1, n (%) | 277 (64.4) | 137 (65.2) | 140 (63.6) |
| Current smoker, n (%) | 177 (41.2) | 83 (39.5) | 94 (42.7) |
| Remission, n (%) | 166 (38.6) | 72 (34.3) | 94 (42.7) |

### Table 2. Key bioinformatics processing and validation results.

| Stage or validation | Result | Reference or denominator |
| --- | --- | --- |
| Supplied genotype matrix | 524,876 unique rsID variants | 430 participants |
| GRCh37 autosomal SNV reconstruction | 492,071 (93.7%) | 32,805 not reconstructed |
| Official GenoPred reference overlap | 147,438 (28.09%) | Supplied rsID variants |
| Frequency-based final orientation | 147,370 | 68 EUR frequency ties excluded |
| PRS-CS reference compatibility audit | 141,038 (95.7%) | Final oriented set; diagnostic only |
| Conservative orientation sensitivity | 953 excluded; 140,085 retained | 141,038 compatible variants; posterior effects fixed |
| Final PLINK2 target | 22 chromosomes; 430 participants | 147,370 variants |
| Genome-wide genotype round trip | 63,369,100 exact comparisons | 0 mismatches |
| Formatted-target allele-aware validation | 0 mismatches | 147,370 variants × 430 participants |
| Bayesian PRS matrix | 8 traits × 430 participants | 0 missing or non-finite values |
| Ancestry outputs | 430 maximum-probability EUR | 418 in strict model-based EUR set |

### Table 3. Primary, robustness and conservative orientation-sensitivity performance comparisons from repeated nested cross-validation.

| Analysis | Participants | Metric | Comparator | Comparator + PRS | Increment | 95% CI lower | 95% CI upper |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Primary Elastic Net | 430 | AUC | 0.6939 | 0.6884 | -0.0056 | -0.0113 | 0.0004 |
| Primary Elastic Net | 430 | Brier score | 0.2113 | 0.2126 | -0.0013 | -0.0031 | 0.0006 |
| Primary Elastic Net | 430 | Log loss | 0.6159 | 0.6183 | -0.0024 | -0.0067 | 0.0022 |
| Strict EUR Elastic Net | 418 | AUC | 0.6934 | 0.6853 | -0.0081 | -0.0130 | -0.0033 |
| Strict EUR Elastic Net | 418 | Brier score | 0.2110 | 0.2135 | -0.0025 | -0.0040 | -0.0011 |
| Strict EUR Elastic Net | 418 | Log loss | 0.6131 | 0.6194 | -0.0063 | -0.0099 | -0.0028 |
| PRS-focused Elastic Net | 430 | AUC | 0.5082 | 0.4925 | -0.0157 | -0.0363 | 0.0051 |
| PRS-focused Elastic Net | 430 | Brier score | 0.2413 | 0.2436 | -0.0023 | -0.0047 | 0.0001 |
| PRS-focused Elastic Net | 430 | Log loss | 0.6770 | 0.6821 | -0.0051 | -0.0102 | 0.0000 |
| Summary-predictor Elastic Net | 430 | AUC | 0.7035 | 0.6935 | -0.0100 | -0.0203 | 0.0002 |
| Summary-predictor Elastic Net | 430 | Brier score | 0.2123 | 0.2145 | -0.0022 | -0.0063 | 0.0019 |
| Summary-predictor Elastic Net | 430 | Log loss | 0.6249 | 0.6286 | -0.0037 | -0.0146 | 0.0076 |
| Random Forest | 430 | AUC | 0.7182 | 0.7112 | -0.0070 | -0.0143 | 0.0002 |
| Random Forest | 430 | Brier score | 0.2069 | 0.2084 | -0.0015 | -0.0031 | 0.0001 |
| Random Forest | 430 | Log loss | 0.6011 | 0.6048 | -0.0037 | -0.0073 | 0.0001 |
| Orientation sensitivity: primary Elastic Net | 430 | AUC | 0.6939 | 0.6894 | -0.0046 | -0.0099 | 0.0011 |
| Orientation sensitivity: summary-predictor Elastic Net | 430 | AUC | 0.7035 | 0.6927 | -0.0109 | -0.0210 | -0.0006 |

### Table 4. Adjusted individual-PRS associations with remission and cross-validated AUC increments.

| Trait | Adjusted OR per SD | OR 95% CI lower | OR 95% CI upper | Wald p-value | BH-FDR p-value | AUC increment | AUC 95% CI lower | AUC 95% CI upper | Positive AUC repeats (of 10) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MDD | 0.935 | 0.767 | 1.139 | 0.503 | 0.607 | -0.0067 | -0.0208 | 0.0073 | 0 |
| ANX | 0.906 | 0.744 | 1.105 | 0.331 | 0.607 | -0.0051 | -0.0255 | 0.0148 | 0 |
| BIP | 0.904 | 0.741 | 1.103 | 0.320 | 0.607 | 0.0015 | -0.0187 | 0.0218 | 5 |
| SCZ | 0.831 | 0.681 | 1.014 | 0.068 | 0.272 | 0.0117 | -0.0224 | 0.0464 | 10 |
| NEUR | 0.939 | 0.772 | 1.143 | 0.531 | 0.607 | -0.0088 | -0.0228 | 0.0047 | 0 |
| INSOM | 0.932 | 0.766 | 1.135 | 0.485 | 0.607 | -0.0003 | -0.0151 | 0.0142 | 5 |
| SWB | 1.011 | 0.828 | 1.233 | 0.916 | 0.916 | -0.0104 | -0.0142 | -0.0068 | 0 |
| EA | 1.206 | 0.987 | 1.474 | 0.067 | 0.272 | 0.0127 | -0.0210 | 0.0459 | 10 |

### Table 5. Predefined PRS–treatment interaction estimates and cross-validated AUC increments.

| Trait | Ratio of odds ratios | Ratio 95% CI lower | Ratio 95% CI upper | Wald p-value | BH-FDR p-value | Group 1 OR | Group 1 OR 95% CI lower | Group 1 OR 95% CI upper | Group 2 OR | Group 2 OR 95% CI lower | Group 2 OR 95% CI upper | AUC increment | Brier improvement | Log-loss improvement | Positive AUC repeats (of 10) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MDD | 1.087 | 0.731 | 1.616 | 0.680 | 0.750 | 0.892 | 0.662 | 1.201 | 0.970 | 0.745 | 1.262 | -0.0122 | -0.0012 | -0.0026 | 0 |
| NEUR | 0.870 | 0.584 | 1.296 | 0.492 | 0.750 | 1.019 | 0.752 | 1.383 | 0.886 | 0.686 | 1.146 | -0.0054 | -0.0010 | -0.0022 | 0 |
| INSOM | 0.937 | 0.630 | 1.395 | 0.750 | 0.750 | 0.966 | 0.721 | 1.293 | 0.905 | 0.692 | 1.184 | -0.0105 | -0.0013 | -0.0029 | 0 |

## Supplementary tables

### Supplementary Table S1. Complete supplied clinical-variable inventory and final predictor disposition.

| Position | Variable | Source section | Expected storage class | Final disposition | Semantic type | Primary model role | Penalty factor | In 46-variable summary set |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | subjectid | Administrative identifiers | text or categorical | Identifier; excluded | Not applicable | excluded | Not applicable | No |
| 2 | Row.names | Administrative identifiers | text or categorical | Identifier; excluded | Not applicable | excluded | Not applicable | No |
| 3 | bloodsampleid.x | Administrative identifiers | text or categorical | Identifier; excluded | Not applicable | excluded | Not applicable | No |
| 4 | drug | Treatment allocation | binary integer | Retained predictor | binary categorical | forced unpenalised | 0 | Yes |
| 5 | educ | Demographics and baseline characteristics | integer | Retained predictor | count numeric | regularised clinical predictor | 1 | Yes |
| 6 | marital.rec | Demographics and baseline characteristics | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 7 | occup.rec | Demographics and baseline characteristics | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 8 | children.rec | Demographics and baseline characteristics | integer | Retained predictor | count numeric | regularised clinical predictor | 1 | Yes |
| 9 | smoker | Demographics and baseline characteristics | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 10 | bmi_wk0 | Demographics and baseline characteristics | continuous numeric | Retained predictor | continuous numeric | regularised clinical predictor | 1 | Yes |
| 11 | age | Demographics and baseline characteristics | integer | Retained predictor | continuous numeric | regularised clinical predictor | 1 | Yes |
| 12 | ageonset | Demographics and baseline characteristics | integer | Retained predictor | continuous numeric | regularised clinical predictor | 1 | Yes |
| 13 | sex | Demographics and baseline characteristics | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 14 | madrs1wk0 | Baseline MADRS items | integer | Retained predictor | ordinal madrs item | regularised clinical predictor | 1 | No |
| 15 | madrs2wk0 | Baseline MADRS items | integer | Retained predictor | ordinal madrs item | regularised clinical predictor | 1 | No |
| 16 | madrs3wk0 | Baseline MADRS items | integer | Retained predictor | ordinal madrs item | regularised clinical predictor | 1 | No |
| 17 | madrs4wk0 | Baseline MADRS items | integer | Retained predictor | ordinal madrs item | regularised clinical predictor | 1 | No |
| 18 | madrs5wk0 | Baseline MADRS items | integer | Retained predictor | ordinal madrs item | regularised clinical predictor | 1 | No |
| 19 | madrs6wk0 | Baseline MADRS items | integer | Retained predictor | ordinal madrs item | regularised clinical predictor | 1 | No |
| 20 | madrs7wk0 | Baseline MADRS items | integer | Retained predictor | ordinal madrs item | regularised clinical predictor | 1 | No |
| 21 | madrs8wk0 | Baseline MADRS items | integer | Retained predictor | ordinal madrs item | regularised clinical predictor | 1 | No |
| 22 | madrs9wk0 | Baseline MADRS items | integer | Retained predictor | ordinal madrs item | regularised clinical predictor | 1 | No |
| 23 | madrs10wk0 | Baseline MADRS items | integer | Retained predictor | ordinal madrs item | regularised clinical predictor | 1 | No |
| 24 | hamd1wk0 | Baseline HAMD_HRSD items | integer | Retained predictor | ordinal hamd item | regularised clinical predictor | 1 | No |
| 25 | hamd2wk0 | Baseline HAMD_HRSD items | integer | Retained predictor | ordinal hamd item | regularised clinical predictor | 1 | No |
| 26 | hamd3wk0 | Baseline HAMD_HRSD items | integer | Retained predictor | ordinal hamd item | regularised clinical predictor | 1 | No |
| 27 | hamd4wk0 | Baseline HAMD_HRSD items | integer | Retained predictor | ordinal hamd item | regularised clinical predictor | 1 | No |
| 28 | hamd5wk0 | Baseline HAMD_HRSD items | integer | Retained predictor | ordinal hamd item | regularised clinical predictor | 1 | No |
| 29 | hamd6wk0 | Baseline HAMD_HRSD items | integer | Retained predictor | ordinal hamd item | regularised clinical predictor | 1 | No |
| 30 | hamd7wk0 | Baseline HAMD_HRSD items | integer | Retained predictor | ordinal hamd item | regularised clinical predictor | 1 | No |
| 31 | hamd8wk0 | Baseline HAMD_HRSD items | integer | Retained predictor | ordinal hamd item | regularised clinical predictor | 1 | No |
| 32 | hamd9wk0 | Baseline HAMD_HRSD items | integer | Retained predictor | ordinal hamd item | regularised clinical predictor | 1 | No |
| 33 | hamd10wk0 | Baseline HAMD_HRSD items | integer | Retained predictor | ordinal hamd item | regularised clinical predictor | 1 | No |
| 34 | hamd11wk0 | Baseline HAMD_HRSD items | integer | Retained predictor | ordinal hamd item | regularised clinical predictor | 1 | No |
| 35 | hamd12wk0 | Baseline HAMD_HRSD items | integer | Retained predictor | ordinal hamd item | regularised clinical predictor | 1 | No |
| 36 | hamd13wk0 | Baseline HAMD_HRSD items | integer | Retained predictor | ordinal hamd item | regularised clinical predictor | 1 | No |
| 37 | hamd14wk0 | Baseline HAMD_HRSD items | integer | Retained predictor | ordinal hamd item | regularised clinical predictor | 1 | No |
| 38 | hamd15wk0 | Baseline HAMD_HRSD items | integer | Retained predictor | ordinal hamd item | regularised clinical predictor | 1 | No |
| 39 | hamd16wk0 | Baseline HAMD_HRSD items | continuous numeric | Retained predictor | ordinal hamd item source fractional | regularised clinical predictor | 1 | No |
| 40 | hamd17wk0 | Baseline HAMD_HRSD items | integer | Retained predictor | ordinal hamd item | regularised clinical predictor | 1 | No |
| 41 | bdi1wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 42 | bdi2wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 43 | bdi3wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 44 | bdi4wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 45 | bdi5wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 46 | bdi6wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 47 | bdi7wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 48 | bdi8wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 49 | bdi9wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 50 | bdi10wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 51 | bdi11wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 52 | bdi12wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 53 | bdi13wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 54 | bdi14wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 55 | bdi15wk0 | Baseline BDI items | integer | Later exact duplicate; excluded | Not applicable | excluded | Not applicable | No |
| 56 | bdi16wk0 | Baseline BDI items | integer | Later exact duplicate; excluded | Not applicable | excluded | Not applicable | No |
| 57 | bdi17wk0 | Baseline BDI items | integer | Later exact duplicate; excluded | Not applicable | excluded | Not applicable | No |
| 58 | bdi18wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 59 | bdi19awk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 60 | bdi19wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 61 | bdi20wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 62 | bdi21wk0 | Baseline BDI items | integer | Retained predictor | ordinal bdi item | regularised clinical predictor | 1 | No |
| 63 | f1score0 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | continuous numeric | regularised clinical predictor | 1 | Yes |
| 64 | f2score0 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | continuous numeric | regularised clinical predictor | 1 | Yes |
| 65 | f3score0 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | continuous numeric | regularised clinical predictor | 1 | Yes |
| 66 | f61score0 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | continuous numeric | regularised clinical predictor | 1 | Yes |
| 67 | f62score0 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | continuous numeric | regularised clinical predictor | 1 | Yes |
| 68 | f63score0 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | continuous numeric | regularised clinical predictor | 1 | Yes |
| 69 | f64score0 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | continuous numeric | regularised clinical predictor | 1 | Yes |
| 70 | f65score0 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | continuous numeric | regularised clinical predictor | 1 | Yes |
| 71 | f66score0 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | continuous numeric | regularised clinical predictor | 1 | Yes |
| 72 | mscore | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | count numeric | regularised clinical predictor | 1 | Yes |
| 73 | melanch | Depression subtypes symptom factors dimensions and SCAN variables | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 74 | atscore | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | count numeric | regularised clinical predictor | 1 | Yes |
| 75 | atyp | Depression subtypes symptom factors dimensions and SCAN variables | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 76 | anxdep | Depression subtypes symptom factors dimensions and SCAN variables | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 77 | anxsomdep | Depression subtypes symptom factors dimensions and SCAN variables | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 78 | k1 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 79 | k2 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 80 | k3 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 81 | k4 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 82 | k5 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 83 | k6 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 84 | k7 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 85 | k8 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 86 | k9 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 87 | k10 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 88 | k11 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 89 | k12 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 90 | k13 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 91 | k14 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 92 | k15 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 93 | k16 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 94 | k17 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 95 | k18 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 96 | k19 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 97 | k20 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 98 | k21 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 99 | k22 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 100 | k23 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 101 | k24 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 102 | k25 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 103 | k26 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 104 | k27 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 105 | k28 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Later exact duplicate; excluded | Not applicable | excluded | Not applicable | No |
| 106 | k29 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Later exact duplicate; excluded | Not applicable | excluded | Not applicable | No |
| 107 | k30 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Retained predictor | ordinal scan item | regularised clinical predictor | 1 | No |
| 108 | k31 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Later exact duplicate; excluded | Not applicable | excluded | Not applicable | No |
| 109 | k32 | Depression subtypes symptom factors dimensions and SCAN variables | integer | Later exact duplicate; excluded | Not applicable | excluded | Not applicable | No |
| 110 | newbleq1wk0 | Stressful life events | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | No |
| 111 | newbleq2wk0 | Stressful life events | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | No |
| 112 | newbleq4wk0 | Stressful life events | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | No |
| 113 | newbleq5wk0 | Stressful life events | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | No |
| 114 | newbleq6wk0 | Stressful life events | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | No |
| 115 | newbleq7wk0 | Stressful life events | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | No |
| 116 | newbleq8wk0 | Stressful life events | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | No |
| 117 | newbleq9wk0 | Stressful life events | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | No |
| 118 | nevents | Stressful life events | integer | Retained predictor | count numeric | regularised clinical predictor | 1 | Yes |
| 119 | eventsbin | Stressful life events | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 120 | mhssri | Antidepressant medication history | integer | Retained predictor | count numeric | regularised clinical predictor | 1 | Yes |
| 121 | everssri | Antidepressant medication history | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 122 | mhtca | Antidepressant medication history | integer | Retained predictor | count numeric | regularised clinical predictor | 1 | Yes |
| 123 | evertca | Antidepressant medication history | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 124 | mhdual | Antidepressant medication history | integer | Retained predictor | count numeric | regularised clinical predictor | 1 | Yes |
| 125 | everdual | Antidepressant medication history | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 126 | mhoth | Antidepressant medication history | integer | Retained predictor | count numeric | regularised clinical predictor | 1 | Yes |
| 127 | everoth | Antidepressant medication history | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 128 | mhantidep | Antidepressant medication history | integer | Retained predictor | count numeric | regularised clinical predictor | 1 | Yes |
| 129 | everantidep | Antidepressant medication history | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 130 | concssri | Antidepressant medication history | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 131 | conctca | Antidepressant medication history | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 132 | concantidep | Antidepressant medication history | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 133 | mcbenzo | Antidepressant medication history | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 134 | mczhypn | Antidepressant medication history | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 135 | mhmirta | Antidepressant medication history | binary integer | Retained predictor | binary categorical | regularised clinical predictor | 1 | Yes |
| 136 | madrs.total | Derived baseline clinical totals | integer | Retained predictor | continuous numeric | regularised clinical predictor | 1 | Yes |
| 137 | hrsd.total | Derived baseline clinical totals | continuous numeric | Retained predictor | continuous numeric | regularised clinical predictor | 1 | Yes |
| 138 | bdi.total | Derived baseline clinical totals | integer | Retained predictor | continuous numeric | regularised clinical predictor | 1 | Yes |
| 139 | hdremit.all | Outcome | binary integer | Outcome; excluded | Not applicable | excluded | Not applicable | No |

### Supplementary Table S2. Source-only participant-identifier and linkage audit.

| Audit item | Result | Interpretation boundary |
| --- | --- | --- |
| Clinical Row.names completeness | 430/430 complete | Source-only identifier audit |
| Clinical Row.names uniqueness | 430/430 unique | Source-only identifier audit |
| bloodsampleid.x completeness | 430/430 complete | Source-only identifier audit |
| bloodsampleid.x uniqueness | 430/430 unique | Source-only identifier audit |
| Row.names and bloodsampleid.x order | Identical for 430/430 rows | Clinical source fields agree |
| Final chromosome PSAM identifiers | 430/430 matched recorded clinical IID mapping | Target files preserve selected mapping |
| Original genotype-to-clinical linkage | Assumed paired 430-row source order | No genotype-side participant identifiers were supplied |

### Supplementary Table S3. Chromosome-level allele-orientation decisions and exact genotype round-trip validation.

| Chromosome | Matched variants | Retained variants | Excluded EUR frequency ties | Retention (%) | Participants | Genotype comparisons | Genotype mismatches | Round-trip pass |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 12,212 | 12,205 | 7 | 99.9 | 430 | 5,248,150 | 0 | Yes |
| 2 | 12,088 | 12,084 | 4 | 100.0 | 430 | 5,196,120 | 0 | Yes |
| 3 | 9,571 | 9,564 | 7 | 99.9 | 430 | 4,112,520 | 0 | Yes |
| 4 | 8,398 | 8,392 | 6 | 99.9 | 430 | 3,608,560 | 0 | Yes |
| 5 | 8,645 | 8,641 | 4 | 100.0 | 430 | 3,715,630 | 0 | Yes |
| 6 | 8,766 | 8,763 | 3 | 100.0 | 430 | 3,768,090 | 0 | Yes |
| 7 | 8,079 | 8,073 | 6 | 99.9 | 430 | 3,471,390 | 0 | Yes |
| 8 | 7,390 | 7,389 | 1 | 100.0 | 430 | 3,177,270 | 0 | Yes |
| 9 | 7,233 | 7,227 | 6 | 99.9 | 430 | 3,107,610 | 0 | Yes |
| 10 | 7,666 | 7,661 | 5 | 99.9 | 430 | 3,294,230 | 0 | Yes |
| 11 | 7,693 | 7,688 | 5 | 99.9 | 430 | 3,305,840 | 0 | Yes |
| 12 | 7,510 | 7,509 | 1 | 100.0 | 430 | 3,228,870 | 0 | Yes |
| 13 | 5,831 | 5,830 | 1 | 100.0 | 430 | 2,506,900 | 0 | Yes |
| 14 | 5,158 | 5,157 | 1 | 100.0 | 430 | 2,217,510 | 0 | Yes |
| 15 | 4,872 | 4,871 | 1 | 100.0 | 430 | 2,094,530 | 0 | Yes |
| 16 | 4,904 | 4,903 | 1 | 100.0 | 430 | 2,108,290 | 0 | Yes |
| 17 | 4,518 | 4,516 | 2 | 100.0 | 430 | 1,941,880 | 0 | Yes |
| 18 | 4,702 | 4,701 | 1 | 100.0 | 430 | 2,021,430 | 0 | Yes |
| 19 | 3,009 | 3,006 | 3 | 99.9 | 430 | 1,292,580 | 0 | Yes |
| 20 | 4,247 | 4,246 | 1 | 100.0 | 430 | 1,825,780 | 0 | Yes |
| 21 | 2,372 | 2,370 | 2 | 99.9 | 430 | 1,019,100 | 0 | Yes |
| 22 | 2,574 | 2,574 | 0 | 100.0 | 430 | 1,106,820 | 0 | Yes |

### Supplementary Table S4. GWAS inputs fixed before outcome modelling, standardisation counts and sample-size strategy.

| Trait | Trait label | Analysis role | Source file | Input rows | Autosomal rows | Sample-size strategy | Outcome used in standardisation |
| --- | --- | --- | --- | --- | --- | --- | --- |
| MDD | Major depressive disorder | primary | pgc-mdd2025_no23andMe_eur_v3-49-24-11.tsv.gz | 7,363,302 | 7,192,042 | Per-variant NEFF supplied by the source | No |
| ANX | Anxiety | exploratory | ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz | 7,236,041 | 7,236,041 | two multiplied by the source neff half field | No |
| BIP | Bipolar disorder | exploratory | bip2024_eur_no23andMe.gz | 6,939,126 | 6,939,126 | two multiplied by the source neff half field | No |
| SCZ | Schizophrenia | exploratory | PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz | 7,659,767 | 7,659,767 | Per-variant NEFF supplied by the source | No |
| NEUR | Neuroticism | primary | sumstats_neuroticism_ctg_format.txt.gz | 10,958,177 | 10,849,319 | Per-variant N supplied by the source | No |
| INSOM | Insomnia | primary | insomnia_ukb2b_EUR_sumstats_20190311_with_chrX_mac_100.txt.gz | 12,663,596 | 12,436,059 | Effective N calculated from 109548 cases and 277440 controls | No |
| SWB | Subjective well-being | exploratory | SWB_Full.txt.gz | 2,268,674 | 2,268,674 | Constant N of 298420 | No |
| EA | Educational attainment | exploratory | EA4_additive_excl_23andMe.txt.gz | 10,985,947 | 10,985,947 | Constant N of 765283; 23andMe excluded | No |

### Supplementary Table S5. Distribution and completeness of the eight final participant-level PRS.

| Trait | N | Minimum | Q1 | Median | Mean | Q3 | Maximum | SD |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MDD | 430 | -2.353 | -0.381 | 0.218 | 0.219 | 0.831 | 2.473 | 0.898 |
| ANX | 430 | -3.000 | -0.607 | 0.054 | 0.039 | 0.715 | 2.778 | 1.000 |
| BIP | 430 | -3.394 | -0.717 | 0.017 | 0.000 | 0.741 | 2.867 | 1.044 |
| SCZ | 430 | -4.054 | -0.590 | 0.110 | 0.064 | 0.770 | 2.491 | 1.055 |
| NEUR | 430 | -2.797 | -0.554 | 0.119 | 0.105 | 0.835 | 2.665 | 1.002 |
| INSOM | 430 | -3.468 | -0.661 | 0.013 | 0.012 | 0.708 | 2.576 | 0.996 |
| SWB | 430 | -2.801 | -0.802 | -0.101 | -0.122 | 0.475 | 3.563 | 0.981 |
| EA | 430 | -3.172 | -1.002 | -0.281 | -0.277 | 0.394 | 3.265 | 1.010 |

### Supplementary Table S6. Pearson correlations among final participant-level PRS.

| Trait | MDD | ANX | BIP | SCZ | NEUR | INSOM | SWB | EA |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MDD | 1 | 0.444 | 0.210 | 0.256 | 0.331 | 0.262 | -0.189 | -0.200 |
| ANX | 0.444 | 1 | 0.203 | 0.358 | 0.327 | 0.040 | -0.103 | -0.091 |
| BIP | 0.210 | 0.203 | 1 | 0.413 | 0.143 | 0.058 | -0.139 | 0.009 |
| SCZ | 0.256 | 0.358 | 0.413 | 1 | 0.198 | 0.024 | -0.175 | -0.039 |
| NEUR | 0.331 | 0.327 | 0.143 | 0.198 | 1 | 0.141 | -0.126 | -0.109 |
| INSOM | 0.262 | 0.040 | 0.058 | 0.024 | 0.141 | 1 | -0.058 | -0.131 |
| SWB | -0.189 | -0.103 | -0.139 | -0.175 | -0.126 | -0.058 | 1 | 0.103 |
| EA | -0.200 | -0.091 | 0.009 | -0.039 | -0.109 | -0.131 | 0.103 | 1 |

### Supplementary Table S7. Ancestry-assignment, model-based keep-file and projected-PC results.

| Output | Population OR reference | Participants | Mean OR number | Range or integrity | Use |
| --- | --- | --- | --- | --- | --- |
| Maximum-probability assignment | EUR | 430 | 0.987 | 0.743–0.999 | Descriptive assignment |
| Model-based keep file | AFR | 0 | – | – | No retained participants |
| Model-based keep file | AMR | 0 | – | – | No retained participants |
| Model-based keep file | CSA | 0 | – | – | No retained participants |
| Model-based keep file | EAS | 0 | – | – | No retained participants |
| Model-based keep file | EUR | 418 | – | – | Strict sensitivity membership |
| Model-based keep file | MID | 0 | – | – | No retained participants |
| Outside model-based keep files | – | 12 | – | – | Retained in primary analysis |
| Projected principal components | 1KG+HGDP trans-ancestry reference | 430 | 6 PCs | Complete and finite | Forced ancestry adjustment covariates |

### Supplementary Table S8. Complete primary, predefined robustness and conservative orientation-sensitivity performance results.

| Analysis | Participants | Metric | Comparator | Comparator + PRS | Increment | 95% CI lower | 95% CI upper |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Primary Elastic Net | 430 | AUC | 0.6939 | 0.6884 | -0.0056 | -0.0113 | 0.0004 |
| Primary Elastic Net | 430 | Brier score | 0.2113 | 0.2126 | -0.0013 | -0.0031 | 0.0006 |
| Primary Elastic Net | 430 | Log loss | 0.6159 | 0.6183 | -0.0024 | -0.0067 | 0.0022 |
| Strict EUR Elastic Net | 418 | AUC | 0.6934 | 0.6853 | -0.0081 | -0.0130 | -0.0033 |
| Strict EUR Elastic Net | 418 | Brier score | 0.2110 | 0.2135 | -0.0025 | -0.0040 | -0.0011 |
| Strict EUR Elastic Net | 418 | Log loss | 0.6131 | 0.6194 | -0.0063 | -0.0099 | -0.0028 |
| PRS-focused Elastic Net | 430 | AUC | 0.5082 | 0.4925 | -0.0157 | -0.0363 | 0.0051 |
| PRS-focused Elastic Net | 430 | Brier score | 0.2413 | 0.2436 | -0.0023 | -0.0047 | 0.0001 |
| PRS-focused Elastic Net | 430 | Log loss | 0.6770 | 0.6821 | -0.0051 | -0.0102 | 0.0000 |
| Summary-predictor Elastic Net | 430 | AUC | 0.7035 | 0.6935 | -0.0100 | -0.0203 | 0.0002 |
| Summary-predictor Elastic Net | 430 | Brier score | 0.2123 | 0.2145 | -0.0022 | -0.0063 | 0.0019 |
| Summary-predictor Elastic Net | 430 | Log loss | 0.6249 | 0.6286 | -0.0037 | -0.0146 | 0.0076 |
| Random Forest | 430 | AUC | 0.7182 | 0.7112 | -0.0070 | -0.0143 | 0.0002 |
| Random Forest | 430 | Brier score | 0.2069 | 0.2084 | -0.0015 | -0.0031 | 0.0001 |
| Random Forest | 430 | Log loss | 0.6011 | 0.6048 | -0.0037 | -0.0073 | 0.0001 |
| Orientation sensitivity: primary Elastic Net | 430 | AUC | 0.6939 | 0.6894 | -0.0046 | -0.0099 | 0.0011 |
| Orientation sensitivity: summary-predictor Elastic Net | 430 | AUC | 0.7035 | 0.6927 | -0.0109 | -0.0210 | -0.0006 |

### Supplementary Table S9. Complete individual-PRS association and predictive-increment results.

| Trait | Adjusted OR per SD | OR 95% CI lower | OR 95% CI upper | Wald p-value | BH-FDR p-value | Positive AUC repeats (of 10) | AUC increment | AUC 95% CI lower | AUC 95% CI upper | Brier increment | Brier lower 95 | Brier upper 95 | Log loss increment | Log loss lower 95 | Log loss upper 95 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MDD | 0.935 | 0.767 | 1.139 | 0.503 | 0.607 | 0 | -0.0067 | -0.0208 | 0.0073 | -0.0011 | -0.0027 | 0.0004 | -0.0023 | -0.0056 | 0.0010 |
| ANX | 0.906 | 0.744 | 1.105 | 0.331 | 0.607 | 0 | -0.0051 | -0.0255 | 0.0148 | -0.0002 | -0.0023 | 0.0020 | -0.0003 | -0.0048 | 0.0043 |
| BIP | 0.904 | 0.741 | 1.103 | 0.320 | 0.607 | 5 | 0.0015 | -0.0187 | 0.0218 | -0.0006 | -0.0028 | 0.0016 | -0.0016 | -0.0063 | 0.0032 |
| SCZ | 0.831 | 0.681 | 1.014 | 0.068 | 0.272 | 10 | 0.0117 | -0.0224 | 0.0464 | 0.0007 | -0.0033 | 0.0048 | 0.0017 | -0.0070 | 0.0105 |
| NEUR | 0.939 | 0.772 | 1.143 | 0.531 | 0.607 | 0 | -0.0088 | -0.0228 | 0.0047 | -0.0008 | -0.0023 | 0.0006 | -0.0018 | -0.0049 | 0.0013 |
| INSOM | 0.932 | 0.766 | 1.135 | 0.485 | 0.607 | 5 | -0.0003 | -0.0151 | 0.0142 | -0.0010 | -0.0027 | 0.0006 | -0.0021 | -0.0058 | 0.0013 |
| SWB | 1.011 | 0.828 | 1.233 | 0.916 | 0.916 | 0 | -0.0104 | -0.0142 | -0.0068 | -0.0011 | -0.0015 | -0.0008 | -0.0024 | -0.0032 | -0.0017 |
| EA | 1.206 | 0.987 | 1.474 | 0.067 | 0.272 | 10 | 0.0127 | -0.0210 | 0.0459 | 0.0003 | -0.0039 | 0.0043 | 0.0008 | -0.0083 | 0.0096 |

### Supplementary Table S10. Complete predefined PRS-by-treatment-group interaction results.

| Trait | Ratio of odds ratios | Ratio 95% CI lower | Ratio 95% CI upper | Wald p-value | BH-FDR p-value | Group 1 OR | Group 1 OR 95% CI lower | Group 1 OR 95% CI upper | Group 2 OR | Group 2 OR 95% CI lower | Group 2 OR 95% CI upper | AUC increment | Brier improvement | Log-loss improvement | Positive AUC repeats (of 10) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MDD | 1.087 | 0.731 | 1.616 | 0.680 | 0.750 | 0.892 | 0.662 | 1.201 | 0.970 | 0.745 | 1.262 | -0.0122 | -0.0012 | -0.0026 | 0 |
| NEUR | 0.870 | 0.584 | 1.296 | 0.492 | 0.750 | 1.019 | 0.752 | 1.383 | 0.886 | 0.686 | 1.146 | -0.0054 | -0.0010 | -0.0022 | 0 |
| INSOM | 0.937 | 0.630 | 1.395 | 0.750 | 0.750 | 0.966 | 0.721 | 1.293 | 0.905 | 0.692 | 1.184 | -0.0105 | -0.0013 | -0.0029 | 0 |

### Supplementary Table S11. Software, resampling and reproducibility record.

| Component | Software | Fixed design or seed | Scope |
| --- | --- | --- | --- |
| GenoPred and PRS-CS scoring | R 4.2.3; Python 3.8.20; GenoPred v2.3.3-11-g75039d7; commit 75039d76b7e1918e5d43c0b3470c94e0c72332b2 | PRS-CS a=1, b=0.5, 1000 iterations, 500 burn-in, thinning=5, beta_std=False, chromosome seed=1, phi=auto | GenoPred scoring runtime |
| GenoPred workflow dependencies | Snakemake 7.32.3; PLINK 1.90b6.21; PLINK2 2.00a5.12LM | Chromosomes 1–22 | Formatting, ancestry, scoring and validation |
| Reconstruction and reference work | R 4.3.0; Python 3.11.6; PLINK2 2.00a5.10 | GRCh37 fixed resources | Variant reconstruction and target development |
| Predictive modelling | R 4.2.3; glmnet 4.1-8; randomForest 4.7-1.1; pROC 1.18.5 | 10 repeats × 10 outer folds; 5 inner folds; seed 20260715 | Nested CV and performance analysis |
| Paired performance uncertainty | R 4.2.3 | 10,000 outcome-stratified participant bootstrap replicates; seed 20260715 | Fixed saved OOF predictions; no refitting |
| Formal calibration and coefficient stability | Python 3.13.5; NumPy 2.3.5; pandas 2.2.3; Numba 0.65.1; statsmodels 0.14.6 | 10,000 paired bootstrap replicates; seed 20260729 | Predictions and 100 primary outer-fit coefficient tables |
| Release reporting | Python >=3.11; pinned packages in pyproject.toml | Release-safe aggregate inputs | Tables, figures and release validation |

### Supplementary Table S12. Trait-specific effect of conservative target-orientation exclusion on PRS weights and participant scores.

| Trait | Primary non-zero weights | Excluded nonzero weights | Pearson correlation | Spearman correlation | Mean absolute difference | Maximum absolute difference | Participants with changed score |
| --- | --- | --- | --- | --- | --- | --- | --- |
| MDD | 123,217 | 828 | 0.998 | 0.998 | 0.146 | 0.342 | 430 |
| ANX | 135,070 | 926 | 0.999 | 0.999 | 0.140 | 0.308 | 430 |
| BIP | 119,348 | 825 | 0.995 | 0.995 | 0.083 | 0.289 | 430 |
| SCZ | 125,395 | 863 | 0.998 | 0.997 | 0.056 | 0.222 | 422 |
| NEUR | 126,229 | 783 | 0.996 | 0.995 | 0.167 | 0.424 | 430 |
| INSOM | 127,825 | 797 | 0.994 | 0.993 | 0.085 | 0.296 | 429 |
| SWB | 131,542 | 885 | 0.998 | 0.998 | 0.046 | 0.170 | 426 |
| EA | 138,758 | 934 | 0.998 | 0.997 | 0.125 | 0.322 | 430 |

### Supplementary Table S13. Trait-specific GWAS and final PRS scoring attrition.

| Trait | Source rows | Autosomal rows | Cleaned harmonised GWAS rows | Target score markers | Zero weight markers | Non-zero posterior weights |
| --- | --- | --- | --- | --- | --- | --- |
| MDD | 7,363,302 | 7,192,042 | 128,669 | 147,370 | 24,153 | 123,217 |
| ANX | 7,236,041 | 7,236,041 | 141,079 | 147,370 | 12,300 | 135,070 |
| BIP | 6,939,126 | 6,939,126 | 124,796 | 147,370 | 28,022 | 119,348 |
| SCZ | 7,659,767 | 7,659,767 | 131,155 | 147,370 | 21,975 | 125,395 |
| NEUR | 10,958,177 | 10,849,319 | 131,957 | 147,370 | 21,141 | 126,229 |
| INSOM | 12,663,596 | 12,436,059 | 133,619 | 147,370 | 19,545 | 127,825 |
| SWB | 2,268,674 | 2,268,674 | 137,094 | 147,370 | 15,828 | 131,542 |
| EA | 10,985,947 | 10,985,947 | 144,996 | 147,370 | 8,612 | 138,758 |

### Supplementary Table S14. Formal repeat-wise calibration of the primary clinical and clinical-plus-PRS models.

| Model or contrast | Calibration measure | Estimate | 95% CI lower | 95% CI upper | Ideal value |
| --- | --- | --- | --- | --- | --- |
| Clinical + ancestry PCs | Calibration-in-the-large | 0.004 | -0.0771 | 0.0873 | 0 |
| Clinical + ancestry PCs + 8 PRS | Calibration-in-the-large | 0.003 | -0.0737 | 0.0823 | 0 |
| Combined minus clinical | Calibration-in-the-large difference | -0.000 | -0.0106 | 0.0095 | 0 |
| Clinical + ancestry PCs | Calibration slope | 0.779 | 0.5582 | 1.0382 | 1 |
| Clinical + ancestry PCs + 8 PRS | Calibration slope | 0.788 | 0.5606 | 1.0493 | 1 |
| Combined minus clinical | Calibration slope difference | 0.009 | -0.0205 | 0.0387 | 0 |

### Supplementary Table S15. Descriptive stability of the eight PRS coefficients across primary combined-model outer fits.

| PRS | Non-zero (all fits) | Non-zero (alpha > 0 fits) | Predominant sign | Sign consistency (%) | Median non-zero coefficient | Q1 non-zero coefficient | Q3 non-zero coefficient |
| --- | --- | --- | --- | --- | --- | --- | --- |
| MDD | 55/100 | 2/47 | Positive | 61.8 | 0.002 | -0.005 | 0.011 |
| ANX | 54/100 | 1/47 | Negative | 51.9 | -0.000 | -0.004 | 0.007 |
| BIP | 56/100 | 3/47 | Negative | 96.4 | -0.016 | -0.022 | -0.009 |
| SCZ | 89/100 | 36/47 | Negative | 100.0 | -0.044 | -0.058 | -0.034 |
| NEUR | 53/100 | 0/47 | Negative | 67.9 | -0.005 | -0.010 | 0.002 |
| INSOM | 55/100 | 2/47 | Positive | 85.5 | 0.012 | 0.003 | 0.025 |
| SWB | 54/100 | 1/47 | Negative | 94.4 | -0.016 | -0.026 | -0.011 |
| EA | 76/100 | 23/47 | Positive | 100.0 | 0.035 | 0.023 | 0.046 |
