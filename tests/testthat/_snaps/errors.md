# an unimplemented estimator is a classed error

    Code
      svysem("y ~ x1", design, estimator = "PML")
    Condition
      Error:
      ! estimator "PML" is not implemented. Available estimators are "MLM", "WLSMV" and "MLR".

# missing = ml names MLR and is refused for MLM

    Code
      svysem("y ~ x1", design, missing = "ml")
    Condition
      Error:
      ! missing = "ml" is not available for MLM. MLM is complete data only. The MLR estimator will handle missing = "ml".

# an object that is not a survey design is refused

    Code
      svysem("y ~ x1", data.frame(y = 1, x1 = 1))
    Condition
      Error:
      ! design must be a survey.design or survey.design2 object. The class supplied was: data.frame.

# a model variable missing from the design is refused

    Code
      svysem("y ~ not_in_design", design)
    Condition
      Error:
      ! Model variable(s) not found in the design: not_in_design.

# negative weights and a design with no positive weight are refused

    Code
      svysem("y ~ x1", negative)
    Condition
      Error:
      ! Sampling weights must be non-negative and finite, and at least one weight must be positive. Zero weights are kept.

# ordered indicators are refused

    Code
      svysem("y1 ~ y2", design, ordered = "y1")
    Condition
      Error:
      ! Ordered indicators are not supported for MLM.

# a group= argument is refused

    Code
      svysem("y ~ x1", design, estimator = "MLR", group = "strat")
    Condition
      Error:
      ! A group= argument is not supported. MLR fits a single group. Stratification comes from the design.

# a cluster= argument is refused

    Code
      svysem("y ~ x1", design, cluster = "psu")
    Condition
      Error:
      ! A cluster= argument is not supported. Clustering comes from the design; passing cluster would double-count primary sampling units.

# multiple groups are refused

    Code
      svysem("y ~ x1", design, ngroups = 2)
    Condition
      Error:
      ! ngroups= is not supported. Name the groups with group=.

# missing other than listwise is refused

    Code
      svysem("y ~ x1", design, missing = "pairwise")
    Condition
      Error:
      ! missing = "pairwise" is not supported for MLM. Use listwise deletion, which is the default. Incomplete rows stay in the design with weight 0 so the primary sampling units are kept.

# conditional.x = TRUE is refused

    Code
      svysem("y ~ x1", design, conditional.x = TRUE, fixed.x = TRUE, meanstructure = TRUE)
    Condition
      Error:
      ! MLM does not support this sample-moment layout. Unrecognized moment(s): y~x1. fixed.x = TRUE is supported when gamma holds the observed second moments, and the means when a mean structure is estimated. conditional.x = TRUE is not supported.

