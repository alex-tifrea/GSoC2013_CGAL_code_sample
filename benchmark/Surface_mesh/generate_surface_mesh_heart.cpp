#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

#include <CGAL/internal/Finite_support_distribution.h>
#include <CGAL/internal/Weighted_random_generator.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

using namespace std;

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Direction_3 Direction_3;
typedef GT::Point_3 Point_3;
typedef GT::Triangle_3 Triangle_3;
typedef GT::Plane_3 Plane_3;
typedef GT::FT FT;

typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

//Random generator
typedef CGAL::Random_points_in_tetrahedron_3<Point_3> PointGen;
typedef CGAL::internal::Weighted_random_generator<PointGen>
	GeneratorWithWeight;

typedef typename C2t3::Vertex_handle Vertex;

FT heart_function (Point_3 p) {  // radius = 2

  return (2*p.x()*p.x()+p.y()*p.y()+p.z()*p.z()-1)*(2*p.x()*p.x()+p.y()*p.y()+p.z()*p.z()-1)*(2*p.x()*p.x()+p.y()*p.y()+p.z()*p.z()-1) - (0.1*p.x()*p.x()+p.y()*p.y())*p.z()*p.z()*p.z();
}

int main()
{
	ofstream f1("heart");
	ofstream f2("points_heart");
	ofstream f3("normals_heart");
	Tr tr;            // 3D-Delaunay triangulation
	C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

	// defining the surface
	Surface_3 surface(heart_function,             // pointer to function
	      	    Sphere_3(CGAL::ORIGIN, 2.)); // bounding sphere
	// Note that "2." above is the *squared* radius of the bounding sphere!

	// defining meshing criteria
	CGAL::Surface_mesh_default_criteria_3<Tr> criteria(30.,  // angular bound
	      					     0.1,  // radius bound
	      					     0.1); // distance bound
	// meshing surface
	CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

	CGAL::Random rand;
	
	int nr = 1000;
	std::vector<Point_3> points;
	points.reserve(nr);
	CGAL::Random_points_on_surface_mesh_3<Point_3, C2t3> g(c2t3);

	CGAL::cpp11::copy_n( g, nr, std::back_inserter(points));

	std::cout << "Actual number of facet: " <<c2t3.number_of_facets() << "\n";
	int cont = 0;
	C2t3::Facet_iterator it = c2t3.facets_begin();
	f1 << c2t3.number_of_facets() << '\n';
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
		f1 << aux.vertex(0).x() << " "<< aux.vertex(0).y() << " "<<
			aux.vertex(0).z() << "\n";
		f1 << aux.vertex(1).x() << " "<< aux.vertex(1).y() << " "<<
			aux.vertex(1).z() << "\n";
		f1 << aux.vertex(2).x() << " "<< aux.vertex(2).y() << " "<<
			aux.vertex(2).z() << "\n";
		cont++;
	}
	std::cout << "Number of facets counted by me " << cont << '\n';

	std::cout << "Normals are: "<<std::endl;
	it = c2t3.facets_begin();
	int count = 0;
	int countV = 0;
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
		Plane_3 pl = aux.supporting_plane();
		Direction_3 N = pl.orthogonal_direction();
//		FT N[3];
//		Point_3 U(aux.vertex(1).x()-aux.vertex(0).x(),
//				aux.vertex(1).y()-aux.vertex(0).y(),
//				aux.vertex(1).z()-aux.vertex(0).z());
//		Point_3 V(aux.vertex(2).x()-aux.vertex(0).x(),
//				aux.vertex(2).y()-aux.vertex(0).y(),
//				aux.vertex(2).z()-aux.vertex(0).z());
//		N[0] = U.y()*V.z() - U.z()*V.y();
//		N[1] = U.z()*V.x() - U.x()*V.z();
//		N[2] = U.x()*V.y() - U.y()*V.x();
//		Point_3 norm(N[0],N[1],N[2]);
//		if (CGAL::orientation(norm, aux.vertex(0), aux.vertex(1),
//					aux.vertex(2)) == CGAL::NEGATIVE) {
//			N[0] = -N[0];
//			N[1] = -N[1];
//			N[2] = -N[2];
//			std::cout << "BANG\n";
//		}
//		f3 << N[0] << " " << N[1] << " " << N[2] << std::endl;

		Point_3 sample(0, 0, 0.25);
		if (CGAL::orientation(sample, aux.vertex(0), aux.vertex(1),
					aux.vertex(2)) == CGAL::NEGATIVE) {
			std::cout << "BANG\n";
			N = -N;
			count++;
		}
		f3 << N.dx() << " " << N.dy() << " " << N.dz() << std::endl;

//		Point_3 sample1(0, 0, 0.25);
//		if (CGAL::orientation(sample1, aux.vertex(0), aux.vertex(1),
//					aux.vertex(2)) == CGAL::NEGATIVE) {
//			std::cout << "BANG\n";
//			countV++;
//		}
	}
//	std::cout << count << " negatives\n";
//	std::cout << countV << " negatives after checking\n";

	std::cout << "The generated points are: " << std::endl;
	f2 << nr << '\n';
	for (int i = 0; i < nr; i++) {
		f2 << points[i].x() << " " << points[i].y() << " " <<
			points[i].z() << std::endl;
	}
	return 0;
}
