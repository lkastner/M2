compute#Fan#computedFVector = method()
compute#Fan#computedFVector Fan := F -> (
   toList apply(0..(dim F), d -> #faces(dim F - d,F))
)

compute#Fan#computedDimension = method()
compute#Fan#computedDimension Fan := F -> (
   R := rays F;
   MC := maxCones F;
   L := linealitySpace F;
   MC = apply(MC, m -> R_m);
   MC = apply(MC, r -> rank (r | L));
   max MC
)


compute#Fan#computedFacesThroughRays = method()
compute#Fan#computedFacesThroughRays Fan := F -> (
   MC := computedMaxCones F;
   raysF := rays F;
   dimF := dim F;
   linealityF := linealitySpace F;
   result := new MutableHashTable;
   for i from 0 to dim F do result#i = {};
   for C in MC do (
      dimC := dim C;
      raysC := rays C;
      facesC := faces C;
      rc := rayCorrespondenceMap(raysC, linealityF, raysF);
      for i in keys facesC do (
         codimInF := i + dimF - dimC;
         codimiCones := facesC#i;
         codimiCones = apply(codimiCones,
            c -> (
               sort apply(c, e -> rc#e)
            )
         );
         result#codimInF = sort unique flatten {result#codimInF, codimiCones};
      );
   );
   return hashTable pairs result
)

-- UNEXPORTED METHOD
-- TODO: export this
-- Returns the smooth cones of a fan as a list of index lists
smoothCones = method()
smoothCones Fan := List => F -> F.cache.smoothCones ??= (
   result := {};
   raysF := rays F;
   linealityF := linealitySpace F;
   cones := getProperty(F, computedFacesThroughRays);
   for i in keys cones do (
      for cone in cones#i do (
         if spanSmoothCone(transpose(raysF_cone), transpose(linealityF)) then (
            result = append(result, cone)
         )
      )
   );
   result
)

compute#Fan#ambientDimension = method()
compute#Fan#ambientDimension Fan := F -> (
   if hasProperty(F, rays) then return numRows rays F
   else if hasProperty(F, computedLinealityBasis) then return numRows linealitySpace F
   else if hasProperty(F, inputRays) then return numRows getProperty(F, inputRays)
   else if hasProperty(F, inputLinealityGenerators) then return numRows getProperty(F, inputLinealityGenerators)
   else error("No property available to compute ambient dimension.")
)

