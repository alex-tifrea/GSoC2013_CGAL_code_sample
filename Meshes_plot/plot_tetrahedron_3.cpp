#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
using namespace std;

int size, points;
float v[4][3];
float **p;
float **n;

void display() {
  glClear(GL_COLOR_BUFFER_BIT);
glEnable(GL_BLEND);
glEnable(GL_LINE_SMOOTH);
glEnable(GL_POINT_SMOOTH);
glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

glPointSize(6);
glEnable(GL_POINT_SPRITE);

glPushMatrix();
glDisable(GL_LIGHTING);
	glColor3f(1,1,0.4);
	glBegin(GL_LINE_LOOP);
	glVertex3f(v[0][0],v[0][1],v[0][2]);
	glVertex3f(v[1][0],v[1][1],v[1][2]);
	glVertex3f(v[2][0],v[2][1],v[2][2]);
	glEnd();

	glColor3f(1,1,0.4);
	glBegin(GL_LINE_LOOP);
	glVertex3f(v[0][0],v[0][1],v[0][2]);
	glVertex3f(v[1][0],v[1][1],v[1][2]);
	glVertex3f(v[3][0],v[3][1],v[3][2]);
	glEnd();

	glColor3f(1,1,0.4);
	glBegin(GL_LINE_LOOP);
	glVertex3f(v[3][0],v[3][1],v[3][2]);
	glVertex3f(v[1][0],v[1][1],v[1][2]);
	glVertex3f(v[2][0],v[2][1],v[2][2]);
	glEnd();

	glColor3f(1,1,0.4);
	glBegin(GL_LINE_LOOP);
	glVertex3f(v[0][0],v[0][1],v[0][2]);
	glVertex3f(v[3][0],v[3][1],v[3][2]);
	glVertex3f(v[2][0],v[2][1],v[2][2]);
	glEnd();


glDisable(GL_LIGHTING);
for (int i=0; i < points; i++){
GLUquadric *qsphere=gluNewQuadric();
gluQuadricOrientation(qsphere,GLU_OUTSIDE);
	glColor3f(1,1,1);
	glPushMatrix();
	glTranslatef(p[i][0], p[i][1], p[i][2]);
	gluSphere(qsphere,0.010,10,10);
	glPopMatrix();
}



glEnable(GL_LIGHTING);
glPopMatrix();

glEnable(GL_LIGHT0);
glEnable(GL_LIGHTING);
float lpos[4] = { -.2f, .2f, .9797958971f, 0.0f };
glLightfv(GL_LIGHT0,GL_POSITION,lpos);
glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0f);

glEnable(GL_NORMALIZE);
glShadeModel(GL_SMOOTH);
glFlush();

}


void init() {
  // Set the current clear color to black and the current drawing color to
  // white.
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glColor3f(1.0, 1.0, 1.0);

  // Set the camera lens to have a 60 degree (vertical) field of view, an
  // aspect ratio of 4/3, and have everything closer than 1 unit to the
  // camera and greater than 40 units distant clipped away.
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0, 4.0/3.0, 1, 90);

  // Position camera at (4, 6, 5) looking at (0, 0, 0) with the vector
  // <0, 1, 0> pointing upward.
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(11, 11, 17, 6, 6, 6, 1, 6, 9);
}

// Initializes GLUT, the display mode, and main window; registers callbacks;
// does application initialization; enters the main event loop.
int main(int argc, char** argv) {
ifstream f1("tetra");
ifstream f2("points_tetra");
f2 >> points;

for (int i = 0; i < 4; i++) {
	f1 >> v[i][0] >> v[i][1] >> v[i][2];
	std::cout <<v[i][0] << " " << v[i][1] << " " << v[i][2] << '\n';
}

p = new float*[points];
for (int i=0;i<points;i++) {
	p[i] = new float[3];
}
for (int i = 0; i < points; i++) {
	f2 >> p[i][0] >> p[i][1] >> p[i][2];
}

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowPosition(80, 80);
  glutInitWindowSize(800, 600);
  glutCreateWindow("A Heart");
  glutDisplayFunc(display);
  init();
  glutMainLoop();
for (int i=0;i<points;i++) {
	delete[] p[i];
}
delete[] p;
  return 0;
}
