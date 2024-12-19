-- PURPOSE : Giving the generating Cones of the Fan
--   INPUT : 'F'  a Fan
--  OUTPUT : a List of Cones
maxCones = method(TypicalValue => List)
maxCones Fan := F -> maxObjects F


--   INPUT : 'F'  a Fan
--  OUTPUT : a Matrix, where the column vectors are a basis of the lineality space
linSpace Fan := F -> linealitySpace F


--   INPUT : 'F',  a Fan
--  OUTPUT : 'true' or 'false'
isPointed Fan := F -> all(values getProperty(F, honestMaxObjects), c->isPointed c)


-- PURPOSE : Checks if the input is smooth
--   INPUT : 'F'  a Fan
--  OUTPUT : 'true' or 'false'
isSmooth Fan := {} >> o -> F -> getProperty(F, smooth)


isPolytopal = method(TypicalValue => Boolean)
isPolytopal Fan := F -> getProperty(F, polytopal)

-- PURPOSE : Giving the k dimensional Cones of the Fan
--   INPUT : (k,F)  where 'k' is a positive integer and F is a Fan 
--  OUTPUT : a List of Cones
cones = method(TypicalValue => List)
cones(ZZ,Fan) := (k,F) -> (
   d := dim F;
   faces := getProperty(F, computedFacesThroughRays);
   faces#(d-k)
)


-- PURPOSE : Returning a polytope of which the fan is the normal if the fan is polytopal
--   INPUT : 'F',  a Fan
--  OUTPUT : A Polytope of which 'F' is the normal fan
polytope = method(TypicalValue => Polyhedron)
polytope Fan := F -> getProperty(F, computedPolytope)



isPure Fan := F -> getProperty(F, pure)
isComplete Fan := F -> getProperty(F, computedComplete)
maxObjects Fan := F -> getProperty(F, generatingObjects)

objectsOfDim(ZZ, Fan) := (k,F) -> (
	-- Checking for input errors
	if k < 0 or dim F < k then error("k must be between 0 and the dimension of the polyhedral object family.");
	L := select(values getProperty(F, honestMaxObjects), C -> dim C >= k);
	-- Collecting the 'k'-dim faces of all generating cones of dimension greater than 'k'
	unique flatten apply(L, C -> faces(dim(C)-k,C)))


isWellDefined Fan := F -> getProperty(F, isWellDefined)


facesAsCones(ZZ, Fan) := (d, F) -> (
   raysF := rays F;
   linF := linealitySpace F;
   result := faces(d, F);
   apply(result, f -> coneFromVData(raysF_f, linF))
)
