///////////////// Definitions of several famous surfaces /////////////////
double sphere_function (double, double, double);  // (c=(0,0,0), r=1)
double ellipsoid_function (double, double, double);  // (c=(0,0,0), r=1)
double torus_function (double, double, double);  // (c=(0,0,0), r=2)
double chair_function (double, double, double);  // (c=(0,0,0), r=6)
double tanglecube_function (double, double, double);  // (c=(0,0,0), r=4)
double octic_function (double, double, double);  // (c=(0,0,0), r=2)
double heart_function (double, double, double);  // (c=(0,0,0), r=2)
double klein_function (double, double, double);  // (c=(0,0,0), r=4)
double ring_function (double, double, double);  // (c=(0,0,0), r=?)
double false_knot_function (double, double, double);  // (c=(0,0,0), r=1)
double knot1_function (double, double, double);  // (c=(0,0,0), r=4)
double knot2_function (double, double, double);  // (c=(0,0,0), r=4)
double knot3_function (double, double, double);  // (c=(0,0,0), r=4)
double cube_function (double, double, double);  // (c=(0,0,0), r=2)


// Sphere (r=1)
template <int Sq_radius>
double sphere_function (double x, double y, double z) // (c=(0,0,0), r=Sq_radius)
{
  double x2=x*x, y2=y*y, z2=z*z;
  return (x2+y2+z2)/Sq_radius - 1;
}

// "Tanglecube" (r=4)
double tanglecube_function (double x, double y, double z) {
  double x2=x*x, y2=y*y, z2=z*z;
  double x4=x2*x2, y4=y2*y2, z4=z2*z2;

  return x4 - 5*x2 + y4 - 5*y2 + z4 - 5*z2 + 11.8;
}

double cube_function (double x, double y, double z){
  if( x < 1 && -x < 1 &&
      y < 1 && -y < 1 &&
      z < 1 && -z < 1 )
    return -1.;
  return 1.;
}

// Barth's octic surface (degree 8)
double octic_function (double x, double y, double z) {  // r=2
  double x2=x*x, y2=y*y, z2=z*z;
  double x4=x2*x2, y4=y2*y2, z4=z2*z2;
  double x6=x4*x2, y6=y4*y2, z6=z4*z2;
  double x8=x4*x4, y8=y4*y4, z8=z4*z4;

  return 43.30495169 *x2*  y2  + 43.30495169 *x2*  z2  + 43.30495169 *y2 *    z2 + 44.36067980 *x6*  y2     + 44.36067980* x6 * z2  + 66.54101970* x4*  y4  + 66.54101970* x4 * z4    + 44.36067980 *x2 * y6  - 11.70820393* x2  - 11.70820393* y2  - 11.70820393* z2     + 37.65247585 *x4 + 37.65247585 *y4  + 37.65247585* z4  + 11.09016995* x8     + 11.09016995  *y8  + 11.09016995* z8  + 133.0820394 *x2*  y4*  z2     + 133.0820394   *x2 * y2 * z4  + 44.36067980* x2 * z6  + 44.36067980 *y6 * z2     + 66.54101970 *y4 * z4  + 44.36067980 *y2 * z6  + 133.0820394*    x4*  y2 * z2     - 91.95742756 *x4 * y2  - 91.95742756 *x4  *z2  - 91.95742756* x2 * y4     - 91.95742756 *x2  *z4  - 91.95742756* y4 *   z2  - 183.9148551* x2  *y2  *z2     - 30.65247585 *x6  - 30.65247585* y6  - 91.95742756 *y2 * z4  - 30.65247585* z6   + 1.618033988;
}

// Torus (r=2)
double torus_function (double x, double y, double z) {
  double x2=x*x, y2=y*y, z2=z*z;
  double x4=x2*x2, y4=y2*y2, z4=z2*z2;

  return x4  + y4  + z4  + 2 *x2*  y2  + 2*
    x2*z2  + 2*y2*  z2  - 5 *x2  + 4* y2  - 5*z2+4;
}

// "Heart"
double heart_function (double x, double y, double z) {  // radius = 2

  return (2*x*x+y*y+z*z-1)*(2*x*x+y*y+z*z-1)*(2*x*x+y*y+z*z-1) - (0.1*x*x+y*y)*z*z*z;
}


// Klein's bottle
double klein_function (double x, double y, double z) {  // radius = 4

  return (x*x+y*y+z*z+2*y-1)*((x*x+y*y+z*z-2*y-1)*(x*x+y*y+z*z-2*y-1)-8*z*z)+16*x*z*(x*x+y*y+z*z-2*y-1);
}

// False knot
double false_knot_function (double x, double y, double z) {  // radius = 1
  double d=1.2, e=0.1;

  double f1 = x*(x*x-z*z)-2*x*z*z-y*y+d*d-x*x-z*z-y*y;

  double m2 = z*(x*x-z*z)+2*x*x*z;
  double f2 = 4*y*y*(d*d-x*x-z*z-y*y) - m2*m2;

  f1 = f1*f1-e*e;
  f2 = f2*f2-e*e;

  if (f1 < 0 && f2 < 0)
    return -1.;
  else if (f1 > 0 || f2 > 0)
    return 1.;
  else
    return 0.;
}


// Knot 1
void puiss(double& x, double& y, int n) {

  double xx = 1, yy = 0;

  while(n>0) {
    if (n&1) {
      double xxx = xx, yyy = yy;
      xx = xxx*x - yyy*y;
      yy = xxx*y + yyy*x;
    }

    double xxx = x, yyy = y;
    x=xxx*xxx-yyy*yyy;
    y=2*xxx*yyy;

    n/=2;
  }

  x = xx;
  y = yy;
}

double knot1_function (double a, double b, double c) {  // radius = 4
  double e=0.1;

  double x, y, z, t, den;
  den=1+a*a+b*b+c*c;
  x=2*a/den;
  y=2*b/den;
  z=2*c/den;
  t=(1-a*a-b*b-c*c)/den;

  double x3=x, y3=y, z2=z, t2=t;
  puiss(x3,y3,3);
  puiss(z2,t2,2);

  double f1 = z2-x3;
  double f2 = t2-y3;

  f1 = f1*f1;
  f2 = f2*f2;
  e=e*e/(den-1);

  if (f1+f2-e < 0)
    return -1.;
  else if (f1+f2-e > 0)
    return 1.;
  else
    return 0.;
}



template <typename FT, typename Point>
class FT_to_point_function_wrapper : public std::unary_function<Point, FT>
{
  typedef FT (*Implicit_function)(FT, FT, FT);
  Implicit_function function;
public:
  FT_to_point_function_wrapper(Implicit_function f) : function(f) {}
  FT operator()(Point p) const { return function(p.x(), p.y(), p.z()); }
};
