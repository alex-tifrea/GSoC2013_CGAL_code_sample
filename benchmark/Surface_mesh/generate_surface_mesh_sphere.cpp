#include <iostream>
#include <CGAL/Timer.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

using namespace std;

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Triangle_3 Triangle_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;

typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

typedef typename C2t3::Vertex_handle Vertex;
	
FT sphere_function (Point_3 p) {
  const FT x2=p.x()*p.x(), y2=p.y()*p.y(), z2=p.z()*p.z();
  return x2+y2+z2-1;
}

int main()
{
	CGAL::Timer t;
	t.start();
	Tr tr;            // 3D-Delaunay triangulation
	C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

	// defining the surface
	Surface_3 surface(sphere_function,             // pointer to function
	      	    Sphere_3(CGAL::ORIGIN, 2.)); // bounding sphere
	// Note that "2." above is the *squared* radius of the bounding sphere!

	// defining meshing criteria
	CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
	      					     0.1,  // radius bound
	      					     0.1); // distance bound
	// meshing surface
	CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

	CGAL::Random rand;
	
	t.stop();
	std::cout << "Time elapsed for building the mesh: " << t.time() << std::endl;
	t.reset();

	t.start();
	int nr = 1000;
	std::vector<Point_3> points;
	points.reserve(nr);
	CGAL::Random_points_on_surface_mesh_3<Point_3, C2t3> g(c2t3);
	t.stop();
	std::cout << "Time elapsed for init Random_points_in_surface_mesh_3: " <<
		t.time() << std::endl;
	t.reset();

	t.start();
	CGAL::cpp11::copy_n( g, nr, std::back_inserter(points));
	t.stop();
	std::cout << "Time elapsed for generating the points: " << t.time() << std::endl;
	t.reset();

//	std::cout << "Coordinates of the vertices are: " <<'\n';
//	Tr::Finite_vertices_iterator iter = tr.finite_vertices_begin();
//	for (; iter != tr.finite_vertices_end(); ++iter) {
//		Point_3 aux = tr.point(iter);
//	//	std::cout << aux.x() << " " << aux.y() << " " << aux.z() << '\n';
//		std::cout << aux.x() << '\n';
//		std::cout << aux.y() << '\n';
//		std::cout << aux.z() << '\n';
//	}

	std::cout << "Actual number of facet: " <<c2t3.number_of_facets() << "\n";
	int cont = 0;
	C2t3::Facet_iterator it = c2t3.facets_begin();
	for (; it != c2t3.facets_end(); ++it) {
		Vertex v[4];
		int k = 0;
		for(int j = 0; j < 4; j++)
		{
			if(j == it->second) continue;
			v[k] = it->first->vertex(j);
			k++;
		}
		Triangle_3 aux(v[0]->point(), v[1]->point(),
				v[2]->point());
		std::cout << aux.vertex(0).x() << " "<< aux.vertex(0).y() << " "<<
			aux.vertex(0).z() << "\n";
		std::cout << aux.vertex(1).x() << " "<< aux.vertex(1).y() << " "<<
			aux.vertex(1).z() << "\n";
		std::cout << aux.vertex(2).x() << " "<< aux.vertex(2).y() << " "<<
			aux.vertex(2).z() << "\n";
		cont++;
	}
	std::cout << "Number of facets counted by me " << cont << '\n';

	std::cout << "The generated points are: " << std::endl;
	for (int i = 0; i < nr; i++) {
		std::cout << points[i].x() << " " << points[i].y() << " " <<
			points[i].z() << std::endl;
	}
	points.clear();
	return 0;
}
