# OmicsBraid: Cross-Omics Effect Inference

## Why OmicsBraid?

Multi-omics studies frequently summarize each layer independently and
then describe RNA/protein/metabolite agreement qualitatively. OmicsBraid
makes that comparison explicit by estimating standardized effects with
uncertainty, accounting for matched-subject dependence, testing joint
multi-omic evidence, quantifying cross-layer heterogeneity, testing
practical equivalence, and evaluating ordered effect trajectories.

Version 0.2.2 retains the hierarchical inference introduced in v0.1.9,
deliberately distinguishes **absence of detected evidence** from
**evidence of practical equivalence**, and separates confirmed braid
patterns from suggestive effect geometry.

## Simulate known biology

``` r

sim <- simulate_braid_data(n_per_group = 60, rho = 0.5, seed = 42)
sim$data
#> <omics_braid_data>
#>  Samples in metadata: 120 
#>  Assays: RNA, Protein, Metabolite 
#>   - RNA: 7 features x 120 samples
#>   - Protein: 7 features x 120 samples
#>   - Metabolite: 7 features x 120 samples
```

## Run the full workflow

For a quick vignette we use fewer bootstrap and Monte Carlo replicates
than are recommended for final analyses. Both the equivalence margin and
the trajectory margin are scientific SESOIs and should be justified in a
real study.

``` r

fit <- run_omics_braid(
  sim$data,
  group = "group",
  reference = "Control",
  comparison = "Disease",
  omic_order = c("RNA", "Protein", "Metabolite"),
  bootstrap_B = 150,
  ci_method = "percentile",
  integrated_ci_method = "percentile",
  ci_min_boot = 75,
  pattern_draws = 500,
  equivalence_margin = 0.30,
  trajectory_margin = 0.15,
  seed = 42
)
fit
#> <omics_braid_result>
#>  Analysis level: feature 
#>  Entities with integrated estimates: 7 
#>  Omic order: RNA -> Protein -> Metabolite 
#>  Pattern counts:
#> 
#>            uncertain        amplification          attenuation 
#>                    2                    1                    1 
#>  concordant_increase            inversion no_detectable_effect 
#>                    1                    1                    1
```

## Layer-specific and integrated inference

``` r

fit$effects
#>           entity       omic     effect        se        var    conf_low
#> 1  concordant_up        RNA  0.9991331 0.1936305 0.03749278  0.68596798
#> 2    attenuation        RNA  1.2550880 0.1997420 0.03989686  0.90352480
#> 3  amplification        RNA  0.4746970 0.1851276 0.03427224  0.08765197
#> 4      inversion        RNA  0.8883029 0.1913666 0.03662117  0.52332491
#> 5      buffering        RNA  1.0137654 0.1939472 0.03761550  0.65890251
#> 6      emergence        RNA -0.1004282 0.1826892 0.03337536 -0.42269332
#> 7           null        RNA  0.1404783 0.1827992 0.03341556 -0.14491433
#> 8  concordant_up    Protein  0.9808467 0.1932406 0.03734192  0.67201246
#> 9    attenuation    Protein  1.0551402 0.1948645 0.03797217  0.68544821
#> 10 amplification    Protein  1.0478368 0.1947002 0.03790817  0.67768971
#> 11     inversion    Protein  0.8905997 0.1914111 0.03663820  0.46540062
#> 12     buffering    Protein  0.0904011 0.1826674 0.03336738 -0.27993766
#> 13     emergence    Protein -0.1346748 0.1827810 0.03340891 -0.50269775
#> 14          null    Protein  0.1863686 0.1829701 0.03347806 -0.09363612
#> 15 concordant_up Metabolite  0.9764470 0.1931477 0.03730604  0.64938039
#> 16   attenuation Metabolite  0.4544151 0.1849154 0.03419372  0.17691150
#> 17 amplification Metabolite  1.1973173 0.1982588 0.03930654  0.87983544
#> 18     inversion Metabolite -1.0484330 0.1947136 0.03791338 -1.40424603
#> 19     buffering Metabolite -0.1070065 0.1827048 0.03338104 -0.47448543
#> 20     emergence Metabolite  1.2203871 0.1988440 0.03953894  0.89688302
#> 21          null Metabolite -0.1235498 0.1827483 0.03339694 -0.46336026
#>     conf_high      p_value n_reference n_comparison  df reference comparison
#> 1   1.3999440 2.469525e-07          60           60 118   Control    Disease
#> 2   1.5996103 3.309351e-10          60           60 118   Control    Disease
#> 3   0.7689709 1.034256e-02          60           60 118   Control    Disease
#> 4   1.2021396 3.452351e-06          60           60 118   Control    Disease
#> 5   1.3911008 1.722654e-07          60           60 118   Control    Disease
#> 6   0.2560018 5.825104e-01          60           60 118   Control    Disease
#> 7   0.4711961 4.421997e-01          60           60 118   Control    Disease
#> 8   1.3409121 3.859092e-07          60           60 118   Control    Disease
#> 9   1.4247137 6.137844e-08          60           60 118   Control    Disease
#> 10  1.3536815 7.374649e-08          60           60 118   Control    Disease
#> 11  1.2207417 3.274390e-06          60           60 118   Control    Disease
#> 12  0.4676582 6.206746e-01          60           60 118   Control    Disease
#> 13  0.2497660 4.612382e-01          60           60 118   Control    Disease
#> 14  0.5713496 3.084051e-01          60           60 118   Control    Disease
#> 15  1.3625680 4.293959e-07          60           60 118   Control    Disease
#> 16  0.8462612 1.399386e-02          60           60 118   Control    Disease
#> 17  1.6492941 1.549140e-09          60           60 118   Control    Disease
#> 18 -0.7545052 7.265111e-08          60           60 118   Control    Disease
#> 19  0.2017524 5.580907e-01          60           60 118   Control    Disease
#> 20  1.6796148 8.387790e-10          60           60 118   Control    Disease
#> 21  0.1791578 4.989992e-01          60           60 118   Control    Disease
#>           p_adj conf_low_analytic conf_high_analytic conf_low_bootstrap
#> 1  5.762226e-07        0.61962425          1.3786419         0.68596798
#> 2  2.316546e-09        0.86360093          1.6465751         0.90352480
#> 3  1.447959e-02        0.11185350          0.8375405         0.08765197
#> 4  6.041614e-06        0.51323124          1.2633745         0.52332491
#> 5  5.762226e-07        0.63363592          1.3938948         0.65890251
#> 6  5.825104e-01       -0.45849254          0.2576361        -0.42269332
#> 7  5.158997e-01       -0.21780164          0.4987582        -0.14491433
#> 8  9.004548e-07        0.60210210          1.3595912         0.67201246
#> 9  2.581127e-07        0.67321282          1.4370676         0.68544821
#> 10 2.581127e-07        0.66623135          1.4294422         0.67768971
#> 11 5.730183e-06        0.51544090          1.2657585         0.46540062
#> 12 6.206746e-01       -0.26762045          0.4484227        -0.27993766
#> 13 5.381112e-01       -0.49291905          0.2235694        -0.50269775
#> 14 4.317671e-01       -0.17224615          0.5449834        -0.09363612
#> 15 7.514428e-07        0.59788441          1.3550095         0.64938039
#> 16 1.959140e-02        0.09198748          0.8168427         0.17691150
#> 17 5.421990e-09        0.80873730          1.5858974         0.87983544
#> 18 1.695193e-07       -1.43006464         -0.6668014        -1.40424603
#> 19 5.580907e-01       -0.46510134          0.2510883        -0.47448543
#> 20 5.421990e-09        0.83066000          1.6101142         0.89688302
#> 21 5.580907e-01       -0.48172984          0.2346303        -0.46336026
#>    conf_high_bootstrap ci_boot_n  ci_method
#> 1            1.3999440       150 percentile
#> 2            1.5996103       150 percentile
#> 3            0.7689709       150 percentile
#> 4            1.2021396       150 percentile
#> 5            1.3911008       150 percentile
#> 6            0.2560018       150 percentile
#> 7            0.4711961       150 percentile
#> 8            1.3409121       150 percentile
#> 9            1.4247137       150 percentile
#> 10           1.3536815       150 percentile
#> 11           1.2207417       150 percentile
#> 12           0.4676582       150 percentile
#> 13           0.2497660       150 percentile
#> 14           0.5713496       150 percentile
#> 15           1.3625680       150 percentile
#> 16           0.8462612       150 percentile
#> 17           1.6492941       150 percentile
#> 18          -0.7545052       150 percentile
#> 19           0.2017524       150 percentile
#> 20           1.6796148       150 percentile
#> 21           0.1791578       150 percentile
fit$integrated
#>          entity n_omics                  omics integrated_effect integrated_se
#> 1 concordant_up       3 RNA;Protein;Metabolite        0.98594971     0.1556271
#> 2   attenuation       3 RNA;Protein;Metabolite        0.86137608     0.1520446
#> 3 amplification       3 RNA;Protein;Metabolite        0.88802608     0.1405536
#> 4     inversion       3 RNA;Protein;Metabolite        0.23862626     0.1492886
#> 5     buffering       3 RNA;Protein;Metabolite        0.19472966     0.1491831
#> 6     emergence       3 RNA;Protein;Metabolite        0.26339346     0.1494749
#> 7          null       3 RNA;Protein;Metabolite        0.07098309     0.1505605
#>       conf_low conf_high   z_value      p_value  W_omnibus df_omnibus
#> 1  0.761730058 1.3096901 6.3353350 2.368261e-10  40.150059          3
#> 2  0.594496249 1.1323254 5.6652839 1.467813e-08  48.784204          3
#> 3  0.646029106 1.1503796 6.3180610 2.648655e-10  51.992915          3
#> 4 -0.064689857 0.4800720 1.5984224 1.099490e-01 108.564869          3
#> 5 -0.121630339 0.4642696 1.3053065 1.917885e-01  46.230380          3
#> 6  0.005799899 0.5779128 1.7621253 7.804813e-02  56.469860          3
#> 7 -0.191007204 0.3479219 0.4714588 6.373132e-01   3.751782          3
#>      p_omnibus      Q_omics df_heterogeneity p_heterogeneity I2_omics
#> 1 9.902563e-09   0.01358846                2    9.932288e-01  0.00000
#> 2 1.450000e-10  16.68876219                2    2.377285e-04 88.01589
#> 3 3.005512e-11  12.07502056                2    2.387496e-03 83.43688
#> 4 2.234444e-23 106.00991457                2    9.555195e-24 98.11338
#> 5 5.066534e-10  44.52655531                2    2.143786e-10 95.50830
#> 6 3.334804e-12  53.36477448                2    2.582180e-12 96.25221
#> 7 2.895447e-01   3.52950911                2    1.712288e-01 43.33490
#>   direction_agreement min_gls_weight max_gls_weight has_negative_gls_weight
#> 1           1.0000000      0.3194338      0.3569288                   FALSE
#> 2           1.0000000      0.2630789      0.4101146                   FALSE
#> 3           1.0000000      0.3154882      0.3611162                   FALSE
#> 4           0.6742767      0.3314425      0.3358423                   FALSE
#> 5           0.6537201      0.1988189      0.4014725                   FALSE
#> 6           0.2968981      0.2850849      0.3715077                   FALSE
#> 7           0.6663352      0.3078795      0.3653995                   FALSE
#>   covariance_condition_number   covariance_mode        p_adj p_omnibus_adj
#> 1                    4.034579 matched_bootstrap 9.270294e-10  1.155299e-08
#> 2                    3.472585 matched_bootstrap 3.424897e-08  2.537499e-10
#> 3                    2.522744 matched_bootstrap 9.270294e-10  7.012861e-11
#> 4                    3.231057 matched_bootstrap 1.539286e-01  1.564111e-22
#> 5                    4.401232 matched_bootstrap 2.237532e-01  7.093148e-10
#> 6                    4.088148 matched_bootstrap 1.365842e-01  1.167181e-11
#> 7                    4.774200 matched_bootstrap 6.373132e-01  2.895447e-01
#>   p_heterogeneity_adj heterogeneity_test_flag has_omnibus_evidence_local
#> 1        9.932288e-01                   FALSE                       TRUE
#> 2        4.160250e-04                    TRUE                       TRUE
#> 3        3.342494e-03                    TRUE                       TRUE
#> 4        6.688636e-23                    TRUE                       TRUE
#> 5        5.002167e-10                    TRUE                       TRUE
#> 6        9.037629e-12                    TRUE                       TRUE
#> 7        1.997669e-01                   FALSE                      FALSE
#>   has_omnibus_evidence evidence_qualified_direction_agreement conf_low_analytic
#> 1                 TRUE                              1.0000000        0.68092623
#> 2                 TRUE                              1.0000000        0.56337405
#> 3                 TRUE                              1.0000000        0.61254613
#> 4                 TRUE                              0.6742767       -0.05397404
#> 5                 TRUE                              0.6537201       -0.09766382
#> 6                 TRUE                              0.2968981       -0.02957192
#> 7                FALSE                                     NA       -0.22411015
#>   conf_high_analytic conf_low_bootstrap conf_high_bootstrap ci_boot_n
#> 1          1.2909732        0.761730058           1.3096901       150
#> 2          1.1593781        0.594496249           1.1323254       150
#> 3          1.1635060        0.646029106           1.1503796       150
#> 4          0.5312266       -0.064689857           0.4800720       150
#> 5          0.4871231       -0.121630339           0.4642696       150
#> 6          0.5563588        0.005799899           0.5779128       150
#> 7          0.3660763       -0.191007204           0.3479219       150
#>    ci_method
#> 1 percentile
#> 2 percentile
#> 3 percentile
#> 4 percentile
#> 5 percentile
#> 6 percentile
#> 7 percentile
fit$equivalence
#>           entity       omic     effect        se        var    conf_low
#> 1  concordant_up        RNA  0.9991331 0.1936305 0.03749278  0.68596798
#> 2    attenuation        RNA  1.2550880 0.1997420 0.03989686  0.90352480
#> 3  amplification        RNA  0.4746970 0.1851276 0.03427224  0.08765197
#> 4      inversion        RNA  0.8883029 0.1913666 0.03662117  0.52332491
#> 5      buffering        RNA  1.0137654 0.1939472 0.03761550  0.65890251
#> 6      emergence        RNA -0.1004282 0.1826892 0.03337536 -0.42269332
#> 7           null        RNA  0.1404783 0.1827992 0.03341556 -0.14491433
#> 8  concordant_up    Protein  0.9808467 0.1932406 0.03734192  0.67201246
#> 9    attenuation    Protein  1.0551402 0.1948645 0.03797217  0.68544821
#> 10 amplification    Protein  1.0478368 0.1947002 0.03790817  0.67768971
#> 11     inversion    Protein  0.8905997 0.1914111 0.03663820  0.46540062
#> 12     buffering    Protein  0.0904011 0.1826674 0.03336738 -0.27993766
#> 13     emergence    Protein -0.1346748 0.1827810 0.03340891 -0.50269775
#> 14          null    Protein  0.1863686 0.1829701 0.03347806 -0.09363612
#> 15 concordant_up Metabolite  0.9764470 0.1931477 0.03730604  0.64938039
#> 16   attenuation Metabolite  0.4544151 0.1849154 0.03419372  0.17691150
#> 17 amplification Metabolite  1.1973173 0.1982588 0.03930654  0.87983544
#> 18     inversion Metabolite -1.0484330 0.1947136 0.03791338 -1.40424603
#> 19     buffering Metabolite -0.1070065 0.1827048 0.03338104 -0.47448543
#> 20     emergence Metabolite  1.2203871 0.1988440 0.03953894  0.89688302
#> 21          null Metabolite -0.1235498 0.1827483 0.03339694 -0.46336026
#>     conf_high      p_value n_reference n_comparison  df reference comparison
#> 1   1.3999440 2.469525e-07          60           60 118   Control    Disease
#> 2   1.5996103 3.309351e-10          60           60 118   Control    Disease
#> 3   0.7689709 1.034256e-02          60           60 118   Control    Disease
#> 4   1.2021396 3.452351e-06          60           60 118   Control    Disease
#> 5   1.3911008 1.722654e-07          60           60 118   Control    Disease
#> 6   0.2560018 5.825104e-01          60           60 118   Control    Disease
#> 7   0.4711961 4.421997e-01          60           60 118   Control    Disease
#> 8   1.3409121 3.859092e-07          60           60 118   Control    Disease
#> 9   1.4247137 6.137844e-08          60           60 118   Control    Disease
#> 10  1.3536815 7.374649e-08          60           60 118   Control    Disease
#> 11  1.2207417 3.274390e-06          60           60 118   Control    Disease
#> 12  0.4676582 6.206746e-01          60           60 118   Control    Disease
#> 13  0.2497660 4.612382e-01          60           60 118   Control    Disease
#> 14  0.5713496 3.084051e-01          60           60 118   Control    Disease
#> 15  1.3625680 4.293959e-07          60           60 118   Control    Disease
#> 16  0.8462612 1.399386e-02          60           60 118   Control    Disease
#> 17  1.6492941 1.549140e-09          60           60 118   Control    Disease
#> 18 -0.7545052 7.265111e-08          60           60 118   Control    Disease
#> 19  0.2017524 5.580907e-01          60           60 118   Control    Disease
#> 20  1.6796148 8.387790e-10          60           60 118   Control    Disease
#> 21  0.1791578 4.989992e-01          60           60 118   Control    Disease
#>           p_adj conf_low_analytic conf_high_analytic conf_low_bootstrap
#> 1  5.762226e-07        0.61962425          1.3786419         0.68596798
#> 2  2.316546e-09        0.86360093          1.6465751         0.90352480
#> 3  1.447959e-02        0.11185350          0.8375405         0.08765197
#> 4  6.041614e-06        0.51323124          1.2633745         0.52332491
#> 5  5.762226e-07        0.63363592          1.3938948         0.65890251
#> 6  5.825104e-01       -0.45849254          0.2576361        -0.42269332
#> 7  5.158997e-01       -0.21780164          0.4987582        -0.14491433
#> 8  9.004548e-07        0.60210210          1.3595912         0.67201246
#> 9  2.581127e-07        0.67321282          1.4370676         0.68544821
#> 10 2.581127e-07        0.66623135          1.4294422         0.67768971
#> 11 5.730183e-06        0.51544090          1.2657585         0.46540062
#> 12 6.206746e-01       -0.26762045          0.4484227        -0.27993766
#> 13 5.381112e-01       -0.49291905          0.2235694        -0.50269775
#> 14 4.317671e-01       -0.17224615          0.5449834        -0.09363612
#> 15 7.514428e-07        0.59788441          1.3550095         0.64938039
#> 16 1.959140e-02        0.09198748          0.8168427         0.17691150
#> 17 5.421990e-09        0.80873730          1.5858974         0.87983544
#> 18 1.695193e-07       -1.43006464         -0.6668014        -1.40424603
#> 19 5.580907e-01       -0.46510134          0.2510883        -0.47448543
#> 20 5.421990e-09        0.83066000          1.6101142         0.89688302
#> 21 5.580907e-01       -0.48172984          0.2346303        -0.46336026
#>    conf_high_bootstrap ci_boot_n  ci_method equiv_margin p_tost_lower
#> 1            1.3999440       150 percentile          0.3 9.775301e-12
#> 2            1.5996103       150 percentile          0.3 3.472336e-15
#> 3            0.7689709       150 percentile          0.3 1.427939e-05
#> 4            1.2021396       150 percentile          0.3 2.656614e-10
#> 5            1.3911008       150 percentile          0.3 6.270789e-12
#> 6            0.2560018       150 percentile          0.3 1.373262e-01
#> 7            0.4711961       150 percentile          0.3 7.984382e-03
#> 8            1.3409121       150 percentile          0.3 1.698458e-11
#> 9            1.4247137       150 percentile          0.3 1.771967e-12
#> 10           1.3536815       150 percentile          0.3 2.216737e-12
#> 11           1.2207417       150 percentile          0.3 2.483885e-10
#> 12           0.4676582       150 percentile          0.3 1.628991e-02
#> 13           0.2497660       150 percentile          0.3 1.828655e-01
#> 14           0.5713496       150 percentile          0.3 3.928120e-03
#> 15           1.3625680       150 percentile          0.3 1.939109e-11
#> 16           0.8462612       150 percentile          0.3 2.253877e-05
#> 17           1.6492941       150 percentile          0.3 2.137554e-14
#> 18          -0.7545052       150 percentile          0.3 9.999394e-01
#> 19           0.2017524       150 percentile          0.3 1.454126e-01
#> 20           1.6796148       150 percentile          0.3 1.035591e-14
#> 21           0.1791578       150 percentile          0.3 1.671380e-01
#>    p_tost_upper    p_tost p_tost_adj p_difference p_difference_adj
#> 1  9.998473e-01 0.9998473  0.9999991 2.469525e-07     5.762226e-07
#> 2  9.999991e-01 0.9999991  0.9999991 3.309351e-10     2.316546e-09
#> 3  8.273275e-01 0.8273275  0.9999991 1.034256e-02     1.447959e-02
#> 4  9.989447e-01 0.9989447  0.9999991 3.452351e-06     6.041614e-06
#> 5  9.998835e-01 0.9998835  0.9999991 1.722654e-07     5.762226e-07
#> 6  1.419501e-02 0.1373262  0.6699840 5.825104e-01     5.825104e-01
#> 7  1.914240e-01 0.1914240  0.6699840 4.421997e-01     5.158997e-01
#> 8  9.997869e-01 0.9997869  0.9999467 3.859092e-07     9.004548e-07
#> 9  9.999467e-01 0.9999467  0.9999467 6.137844e-08     2.581127e-07
#> 10 9.999387e-01 0.9999387  0.9999467 7.374649e-08     2.581127e-07
#> 11 9.989840e-01 0.9989840  0.9999467 3.274390e-06     5.730183e-06
#> 12 1.256010e-01 0.1256010  0.6236705 6.206746e-01     6.206746e-01
#> 13 8.700635e-03 0.1828655  0.6236705 4.612382e-01     5.381112e-01
#> 14 2.672873e-01 0.2672873  0.6236705 3.084051e-01     4.317671e-01
#> 15 9.997693e-01 0.9997693  0.9999982 4.293959e-07     7.514428e-07
#> 16 7.981574e-01 0.7981574  0.9999982 1.399386e-02     1.959140e-02
#> 17 9.999970e-01 0.9999970  0.9999982 1.549140e-09     5.421990e-09
#> 18 2.176606e-12 0.9999394  0.9999982 7.265111e-08     1.695193e-07
#> 19 1.295117e-02 0.1454126  0.5849829 5.580907e-01     5.580907e-01
#> 20 9.999982e-01 0.9999982  0.9999982 8.387790e-10     5.421990e-09
#> 21 1.023370e-02 0.1671380  0.5849829 4.989992e-01     5.580907e-01
#>    equivalent_local equivalent_adjusted state_local state_adjusted equivalent
#> 1             FALSE               FALSE    positive       positive      FALSE
#> 2             FALSE               FALSE    positive       positive      FALSE
#> 3             FALSE               FALSE    positive       positive      FALSE
#> 4             FALSE               FALSE    positive       positive      FALSE
#> 5             FALSE               FALSE    positive       positive      FALSE
#> 6             FALSE               FALSE   uncertain      uncertain      FALSE
#> 7             FALSE               FALSE   uncertain      uncertain      FALSE
#> 8             FALSE               FALSE    positive       positive      FALSE
#> 9             FALSE               FALSE    positive       positive      FALSE
#> 10            FALSE               FALSE    positive       positive      FALSE
#> 11            FALSE               FALSE    positive       positive      FALSE
#> 12            FALSE               FALSE   uncertain      uncertain      FALSE
#> 13            FALSE               FALSE   uncertain      uncertain      FALSE
#> 14            FALSE               FALSE   uncertain      uncertain      FALSE
#> 15            FALSE               FALSE    positive       positive      FALSE
#> 16            FALSE               FALSE    positive       positive      FALSE
#> 17            FALSE               FALSE    positive       positive      FALSE
#> 18            FALSE               FALSE    negative       negative      FALSE
#> 19            FALSE               FALSE   uncertain      uncertain      FALSE
#> 20            FALSE               FALSE    positive       positive      FALSE
#> 21            FALSE               FALSE   uncertain      uncertain      FALSE
#>        state state_basis
#> 1   positive       local
#> 2   positive       local
#> 3   positive       local
#> 4   positive       local
#> 5   positive       local
#> 6  uncertain       local
#> 7  uncertain       local
#> 8   positive       local
#> 9   positive       local
#> 10  positive       local
#> 11  positive       local
#> 12 uncertain       local
#> 13 uncertain       local
#> 14 uncertain       local
#> 15  positive       local
#> 16  positive       local
#> 17  positive       local
#> 18  negative       local
#> 19 uncertain       local
#> 20  positive       local
#> 21 uncertain       local
fit$trend
#>          entity n_layers_trend            trend_omics aligned_direction
#> 1 concordant_up              3 RNA;Protein;Metabolite                 1
#> 2   attenuation              3 RNA;Protein;Metabolite                 1
#> 3 amplification              3 RNA;Protein;Metabolite                 1
#> 4     inversion              3 RNA;Protein;Metabolite                 1
#> 5     buffering              3 RNA;Protein;Metabolite                 1
#> 6     emergence              3 RNA;Protein;Metabolite                 1
#> 7          null              3 RNA;Protein;Metabolite                 1
#>   same_sign_estimates trend_intercept trend_slope trend_slope_se trend_conf_low
#> 1                TRUE       0.9966243 -0.01104217     0.10131680     -0.2096195
#> 2                TRUE       1.3186526 -0.39865933     0.10145449     -0.5975065
#> 3                TRUE       0.5326097  0.37240866     0.11254276      0.1518289
#> 4               FALSE       1.1771754 -0.93562333     0.10686338     -1.1450717
#> 5               FALSE       0.8945462 -0.58189365     0.09207753     -0.7623623
#> 6               FALSE      -0.2993101  0.59755433     0.09907455      0.4033718
#> 7               FALSE       0.2030749 -0.12964903     0.08458881     -0.2954401
#>   trend_conf_high p_trend_zero p_meaningful_amplification
#> 1      0.18753511 9.132131e-01               9.440252e-01
#> 2     -0.19981217 8.514399e-05               1.000000e+00
#> 3      0.59298841 9.361596e-04               2.406523e-02
#> 4     -0.72617495 2.035213e-18               1.000000e+00
#> 5     -0.40142501 2.622319e-10               1.000000e+00
#> 6      0.79173688 1.625850e-09               3.130931e-06
#> 7      0.03614199 1.253505e-01               9.995268e-01
#>   p_meaningful_attenuation p_flat_tost trajectory_margin covariance_mode_trend
#> 1             9.148932e-01  0.08510676              0.15     matched_bootstrap
#> 2             7.124097e-03  0.99287590              0.15     matched_bootstrap
#> 3             9.999983e-01  0.97593477              0.15     matched_bootstrap
#> 4             9.787955e-14  1.00000000              0.15     matched_bootstrap
#> 5             1.362401e-06  0.99999864              0.15     matched_bootstrap
#> 6             1.000000e+00  0.99999687              0.15     matched_bootstrap
#> 7             5.950624e-01  0.40493763              0.15     matched_bootstrap
#>   trend_state_local p_meaningful_amplification_adj p_meaningful_attenuation_adj
#> 1         uncertain                   1.000000e+00                 1.000000e+00
#> 2       attenuation                   1.000000e+00                 1.662289e-02
#> 3     amplification                   8.422831e-02                 1.000000e+00
#> 4    not_applicable                   1.000000e+00                 6.851568e-13
#> 5    not_applicable                   1.000000e+00                 4.768405e-06
#> 6    not_applicable                   2.191651e-05                 1.000000e+00
#> 7    not_applicable                   1.000000e+00                 1.000000e+00
#>   p_flat_tost_adj trend_state_adjusted
#> 1       0.5957473            uncertain
#> 2       1.0000000          attenuation
#> 3       1.0000000            uncertain
#> 4       1.0000000       not_applicable
#> 5       1.0000000       not_applicable
#> 6       1.0000000       not_applicable
#> 7       1.0000000       not_applicable
fit$classification
#>          entity              pattern      pattern_status  suggestive_pattern
#> 1 concordant_up  concordant_increase direction_confirmed concordant_increase
#> 2   attenuation          attenuation           confirmed         attenuation
#> 3 amplification        amplification           confirmed       amplification
#> 4     inversion            inversion           confirmed           inversion
#> 5     buffering            uncertain          unresolved           buffering
#> 6     emergence            uncertain          unresolved           emergence
#> 7          null no_detectable_effect         no_evidence     null_equivalent
#>                                                         interpretation_label
#> 1 concordant_increase (trajectory unresolved; geometry: concordant_increase)
#> 2                                                                attenuation
#> 3                                                              amplification
#> 4                                                                  inversion
#> 5                                                       suggestive_buffering
#> 6                                                       suggestive_emergence
#> 7                           no_detectable_effect (geometry: null_equivalent)
#>                                                              confirmatory_basis
#> 1 all observed layers support one direction; magnitude trajectory is unresolved
#> 2                 same-direction layers plus meaningful negative GLS trajectory
#> 3                 same-direction layers plus meaningful positive GLS trajectory
#> 4                  statistically supported effects occur in opposite directions
#> 5 partial or transient layer pattern does not satisfy a confirmatory braid rule
#> 6 partial or transient layer pattern does not satisfy a confirmatory braid rule
#> 7             joint-null omnibus test not rejected; equivalence not established
#>   n_layers_observed n_positive n_negative n_equivalent n_uncertain
#> 1                 3          3          0            0           0
#> 2                 3          3          0            0           0
#> 3                 3          3          0            0           0
#> 4                 3          2          1            0           0
#> 5                 3          1          0            0           2
#> 6                 3          1          0            0           2
#> 7                 3          0          0            0           3
#>   omnibus_p_local omnibus_evidence_local trend_slope trend_slope_se
#> 1    9.902563e-09                   TRUE -0.01104217     0.10131680
#> 2    1.450000e-10                   TRUE -0.39865933     0.10145449
#> 3    3.005512e-11                   TRUE  0.37240866     0.11254276
#> 4    2.234444e-23                   TRUE -0.93562333     0.10686338
#> 5    5.066534e-10                   TRUE -0.58189365     0.09207753
#> 6    3.334804e-12                   TRUE  0.59755433     0.09907455
#> 7    2.895447e-01                  FALSE -0.12964903     0.08458881
#>   trend_conf_low trend_conf_high p_trend_zero p_meaningful_amplification
#> 1     -0.2096195      0.18753511 9.132131e-01               9.440252e-01
#> 2     -0.5975065     -0.19981217 8.514399e-05               1.000000e+00
#> 3      0.1518289      0.59298841 9.361596e-04               2.406523e-02
#> 4     -1.1450717     -0.72617495 2.035213e-18               1.000000e+00
#> 5     -0.7623623     -0.40142501 2.622319e-10               1.000000e+00
#> 6      0.4033718      0.79173688 1.625850e-09               3.130931e-06
#> 7     -0.2954401      0.03614199 1.253505e-01               9.995268e-01
#>   p_meaningful_attenuation p_flat_tost trajectory_margin trend_state_local
#> 1             9.148932e-01  0.08510676              0.15         uncertain
#> 2             7.124097e-03  0.99287590              0.15       attenuation
#> 3             9.999983e-01  0.97593477              0.15     amplification
#> 4             9.787955e-14  1.00000000              0.15    not_applicable
#> 5             1.362401e-06  0.99999864              0.15    not_applicable
#> 6             1.000000e+00  0.99999687              0.15    not_applicable
#> 7             5.950624e-01  0.40493763              0.15    not_applicable
#>   p_meaningful_amplification_adj p_meaningful_attenuation_adj p_flat_tost_adj
#> 1                   1.000000e+00                 1.000000e+00       0.5957473
#> 2                   1.000000e+00                 1.662289e-02       1.0000000
#> 3                   8.422831e-02                 1.000000e+00       1.0000000
#> 4                   1.000000e+00                 6.851568e-13       1.0000000
#> 5                   1.000000e+00                 4.768405e-06       1.0000000
#> 6                   2.191651e-05                 1.000000e+00       1.0000000
#> 7                   1.000000e+00                 1.000000e+00       1.0000000
#>   trend_state_adjusted covariance_mode_trend uncertainty_mode_pattern
#> 1            uncertain     matched_bootstrap      concordant_increase
#> 2          attenuation     matched_bootstrap              attenuation
#> 3            uncertain     matched_bootstrap            amplification
#> 4       not_applicable     matched_bootstrap                inversion
#> 5       not_applicable     matched_bootstrap                buffering
#> 6       not_applicable     matched_bootstrap                emergence
#> 7       not_applicable     matched_bootstrap          null_equivalent
#>   uncertainty_mode_support deterministic_pattern_support
#> 1                    0.852                         0.852
#> 2                    0.790                         0.790
#> 3                    0.820                         0.820
#> 4                    1.000                         1.000
#> 5                    0.828                            NA
#> 6                    0.718                            NA
#> 7                    0.488                            NA
#>   suggestive_pattern_support pattern_entropy mode_agrees_with_deterministic
#> 1                      0.852       0.4635762                           TRUE
#> 2                      0.790       0.7414827                           TRUE
#> 3                      0.820       0.4899206                           TRUE
#> 4                      1.000       0.0000000                           TRUE
#> 5                      0.828       0.4703161                             NA
#> 6                      0.718       0.6203956                             NA
#> 7                      0.488       0.6234601                             NA
#>   classification_stability
#> 1           direction_high
#> 2                 moderate
#> 3                     high
#> 4                     high
#> 5               unresolved
#> 6               unresolved
#> 7             geometry_low
```

The integrated table separates the multivariate all-zero omnibus test
from the GLS consensus effect. This is important when opposite layer
effects cancel in the consensus.

## Robust interval reporting

Version 0.2.2 can display matched subject-bootstrap intervals without
changing the analytic SE/p-value inferential machinery. `fit$effects`
and `fit$integrated` retain analytic interval columns alongside the
selected percentile bootstrap intervals, while `fit$effect_intervals`
and `fit$consensus_intervals` contain bootstrap diagnostics. BCa layer
intervals are available for targeted/prefiltered analyses through
`bootstrap_effect_intervals(..., method = "bca")`.

``` r

head(fit$effect_intervals)
#>          entity omic  ci_method conf_level boot_n   boot_mean   boot_sd
#> 1 concordant_up  RNA percentile       0.95    150  1.04352743 0.1900172
#> 2   attenuation  RNA percentile       0.95    150  1.26530512 0.1727068
#> 3 amplification  RNA percentile       0.95    150  0.48455419 0.1760240
#> 4     inversion  RNA percentile       0.95    150  0.84351771 0.1840193
#> 5     buffering  RNA percentile       0.95    150  1.01251424 0.1942016
#> 6     emergence  RNA percentile       0.95    150 -0.09846206 0.1896144
#>      boot_bias conf_low_boot conf_high_boot bca_z0 bca_acceleration
#> 1  0.044394338    0.68596798      1.3999440     NA               NA
#> 2  0.010217105    0.90352480      1.5996103     NA               NA
#> 3  0.009857210    0.08765197      0.7689709     NA               NA
#> 4 -0.044785162    0.52332491      1.2021396     NA               NA
#> 5 -0.001251129    0.65890251      1.3911008     NA               NA
#> 6  0.001966149   -0.42269332      0.2560018     NA               NA
#>   bca_prob_low bca_prob_high
#> 1        0.025         0.975
#> 2        0.025         0.975
#> 3        0.025         0.975
#> 4        0.025         0.975
#> 5        0.025         0.975
#> 6        0.025         0.975
head(fit$consensus_intervals)
#>          entity  ci_method conf_level boot_n boot_mean   boot_sd    boot_bias
#> 1 concordant_up percentile       0.95    150 1.0225248 0.1526545  0.036575103
#> 2   attenuation percentile       0.95    150 0.8766985 0.1440815  0.015322377
#> 3 amplification percentile       0.95    150 0.9074674 0.1315093  0.019441312
#> 4     inversion percentile       0.95    150 0.2101759 0.1450974 -0.028450368
#> 5     buffering percentile       0.95    150 0.2010555 0.1556106  0.006325826
#> 6     emergence percentile       0.95    150 0.2663238 0.1561983  0.002930297
#>   conf_low_boot conf_high_boot
#> 1   0.761730058      1.3096901
#> 2   0.594496249      1.1323254
#> 3   0.646029106      1.1503796
#> 4  -0.064689857      0.4800720
#> 5  -0.121630339      0.4642696
#> 6   0.005799899      0.5779128
```

## Ordered GLS trend test

`fit$trend` reports a covariance-aware slope over the prespecified omic
order. For same-direction layer effects, OmicsBraid only calls
attenuation or amplification when the slope exceeds the user-specified
practical trajectory margin with inferential support. Otherwise the
broader concordant direction can be retained with
`pattern_status = "direction_confirmed"`.

## Null-like outcomes

`null_equivalent` means all observed layers passed practical-equivalence
testing. `no_detectable_effect` means the joint-null omnibus test was
not rejected, but the data were not precise enough to establish
equivalence. These are different scientific conclusions.

## Evidence Forest

``` r

plot_evidence_forest(fit, "inversion")
```

![](OmicsBraid-introduction_files/figure-html/forest-1.png)

## Effect Braid

``` r

plot_effect_braid(fit, "inversion")
```

![](OmicsBraid-introduction_files/figure-html/braid-1.png)

## Interpretation

A significant `Q_omics` indicates that the effect vector is poorly
summarized by a single common effect after accounting for estimated
sampling covariance. The I2-like index is descriptive and has no
universal low/moderate/high cutoffs in this development version.

Braid order is descriptive/inferential, not causal. A pattern such as
`inversion` describes supported effects across a prespecified layer
order; it does not prove molecular propagation. Likewise,
`suggestive_pattern` describes effect geometry when confirmatory
precision is insufficient and must not be reported as a confirmed
result.
