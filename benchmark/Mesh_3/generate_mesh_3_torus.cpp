#include <iostream>
#include <CGAL/Timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Mesh_3/implicit_to_labeled_function_wrapper.h>
#include <CGAL/Mesh_3/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include "implicit_functions.h"

#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

using namespace std;
using namespace CGAL::parameters;

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef FT_to_point_function_wrapper<K::FT, K::Point_3> Function;
typedef CGAL::Mesh_3::Implicit_vector_to_labeled_function_wrapper<Function, K>
                                                        Function_wrapper;
typedef Function_wrapper::Function_vector Function_vector;
typedef CGAL::Mesh_3::Labeled_mesh_domain_3<Function_wrapper, K> Mesh_domain;

typedef K::Point_3 Point;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Mesh Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
typedef Mesh_criteria::Facet_criteria    Facet_criteria;
typedef Mesh_criteria::Cell_criteria     Cell_criteria;

typedef Tr::Geom_traits GT;
typedef GT::Tetrahedron_3 Tetrahedron3;

typedef CGAL::FasterMemoryExpensiveTag FastPolicy;
typedef CGAL::SlowerMemoryEfficientTag SlowPolicy;

int main()
{
	CGAL::Timer t;
	t.start();
	// Define functions
	Function f1(&torus_function);
	Function f2(&sphere_function<5>);
	Function f3(&tanglecube_function);
	Function f4(&heart_function);
	Function f5(&klein_function);
	Function f6(&false_knot_function);
	Function f7(&knot1_function);
	Function f8(&octic_function);

	Function_vector v;
	v.push_back(&f1);
	//v.push_back(&f2);
	//v.push_back(&f3);
	//v.push_back(&f4);
	//v.push_back(&f5);
	//v.push_back(&f6);
	//v.push_back(&f7);
	//v.push_back(&f8);

	// Domain (Warning: Sphere_3 constructor uses square radius !)
	Mesh_domain domain(v, K::Sphere_3(CGAL::ORIGIN, 5.*5.), 1e-6);

	// Set mesh criteria
	Facet_criteria facet_criteria(30, 0.2, 0.02); // angle, size, approximation
	Cell_criteria cell_criteria(2., 0.4); // radius-edge ratio, size
	Mesh_criteria criteria(facet_criteria, cell_criteria);

	// Mesh generation
	C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude(), no_perturb());

	// Perturbation (maximum cpu time: 10s, targeted dihedral angle: default)
	CGAL::perturb_mesh_3(c3t3, domain, time_limit = 10);
	
	// Exudation
	CGAL::exude_mesh_3(c3t3,12);

	CGAL::Random rand;

	t.stop();
	std::cout << "Time elapsed for building the mesh: " << t.time() << std::endl;
	t.reset();
	
	t.start();
	int nr = 1500;
	std::vector<Point> points;
	points.reserve(nr);
	CGAL::Random_points_in_mesh_3<Point, C3t3, FastPolicy> g(c3t3);
	t.stop();
	std::cout << "Time elapsed for init Random_points_in_mesh_3: " <<
		t.time() << std::endl;
	t.reset();

	t.start();
	CGAL::cpp11::copy_n( g, nr, std::back_inserter(points));
	t.stop();
	std::cout << "Time elapsed for generating the points: " << t.time() << std::endl;
	t.reset();

	Tr tr = c3t3.triangulation();
	Tetrahedron3 tet;
	Tr::Finite_cells_iterator it = tr.finite_cells_begin();
	for (; it != tr.finite_cells_end(); it++) {
		if (c3t3.is_in_complex(it)) {
			tet = tr.tetrahedron(it);
			break;
		}
	}

	t.start();
	std::vector<Point> points_tet;
	points_tet.reserve(nr);
	CGAL::Random_points_in_tetrahedron_3<Point> g1(tet);
	CGAL::cpp11::copy_n(g1, nr, std::back_inserter(points_tet));
	t.stop();
	std::cout << "Time elapsed for " << nr << " points in one tetrahedron: " <<
		t.time() << std::endl;
	t.reset();

	std::cout << "The generated points are: " << std::endl;
	for (int i = 0; i < nr; i++) {
		std::cout << points[i].x() << " " << points[i].y() << " " <<
			points[i].z() << std::endl;
	}

	points_tet.clear();
	points.clear();
	v.clear();
	return 0;
}

