#include <polymake/Main.h>
#include <polymake/Matrix.h>
#include <polymake/SparseMatrix.h>
#include <polymake/Rational.h>
using namespace polymake;

int rawPolymakeConvexHull() {
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
  return 0; 
}


