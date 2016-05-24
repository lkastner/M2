-----------------------------------------------------------------------------
-- polymake interface
-----------------------------------------------------------------------------


export rawPolymakeConvexHull(e:Expr):Expr := (
     when e is M:RawMutableMatrixCell do possibleEngineError( Ccode(bool, "IM2_SmithNormalForm(", M.p, ")" ) )
     else WrongArgMutableMatrix());
setupfun("rawPolymakeConvexHull",rawPolymakeConvexHull);


