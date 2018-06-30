*** This is a test input file for asfem

[mesh]
  type=asfem
  dim=2
  xmin=0.0
  xmax=1.0
  ymin=0.0
  ymax=1.0
  nx=5
  ny=5
  meshtype=quad4
[]

[variables]
name=phi
[]

[kernels]
  [ex1]
    name=poisson
    variable=phi
    mate=poisson
  []
[]

[materials]
  [poisson]
    variable=phi
  []
[]

