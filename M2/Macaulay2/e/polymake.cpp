#include "engine.h"

#include "relem.hpp"
#include "ring.hpp"
#include "GF.hpp"

#include "coeffrings.hpp"
#include "aring-zz-gmp.hpp"
#include "mat.hpp"
#include "fractionfreeLU.hpp"
#include "LLL.hpp"
#include "exceptions.hpp"

#include "matrix.hpp"
#include "aring-zzp-ffpack.hpp"
#include "mutablemat.hpp"

#include "finalize.hpp"


#include <polymake/Main.h>
#include <polymake/Matrix.h>
#include <polymake/SparseMatrix.h>
#include <polymake/Rational.h>
using namespace polymake;

M2_bool x_rawPolymakeConvexHull(MutableMatrix *M) {
  try {
    const int dim = 4;
    Main pm;
    pm.set_application("polytope");
    perl::Object p("Polytope<Rational>");
    p.take("VERTICES") << (ones_vector<Rational>() | 
       3*unit_matrix<Rational>(dim));
    const Matrix<Rational> f = p.give("FACETS");
    const Vector<Integer> h = p.give("H_STAR_VECTOR");
    cout << "facets" << endl << f << endl << "h* " << h << endl;
  } catch (const std::exception& ex) {
    std::cerr << "ERROR: " << ex.what() << endl; return 1;
  }
  return true; 
}


