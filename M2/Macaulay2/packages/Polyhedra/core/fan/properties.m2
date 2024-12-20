


compute#Fan#ambientDimension = method()
compute#Fan#ambientDimension Fan := F -> (
   if hasProperty(F, rays) then return numRows rays F
   else if hasProperty(F, computedLinealityBasis) then return numRows linealitySpace F
   else if hasProperty(F, inputRays) then return numRows getProperty(F, inputRays)
   else if hasProperty(F, inputLinealityGenerators) then return numRows getProperty(F, inputLinealityGenerators)
   else error("No property available to compute ambient dimension.")
)

