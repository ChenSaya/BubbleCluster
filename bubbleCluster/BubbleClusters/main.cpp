/*
 *  main.c
 *  This program draws several overlapping filled polygons
 *  to demonstrate the effect order has on alpha blending results.
 *  Use the 't' key to toggle the order of drawing polygons.
 */

#include <stdlib.h>
#include <cstdlib>
#include <GL/glut.h>
#include "Bubble.h"
#include <math.h>
	

BubbleGameEngine bubbleGame;
int x_before, y_before;

/*  Initialize alpha blending function.
 */
static void init(void)
{
   glEnable (GL_BLEND);
   glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
   glShadeModel (GL_FLAT);
   glClearColor (1.0, 1.0, 1.0, 0.0);
   
   // game initialize

   bubbleGame.init();
   //bubbleGame.clusteringInitial();

   //bubbleGame.clustering_spliting();
   //bubbleGame.boundingBox();
   //bubbleGame.drawCurve();
}


void display(void)
{
   glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

   bubbleGame.draw();

   glutSwapBuffers();
}

void reshape(int w, int h)
{
   glViewport(0, 0, (GLsizei) w, (GLsizei) h);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();

   gluOrtho2D (0.0, w, 0.0, h);

   /*if (w <= h) 
      gluOrtho2D (0.0, 1.0, 0.0, 1.0*(GLfloat)h/(GLfloat)w);
   else 
      gluOrtho2D (0.0, 1.0*(GLfloat)w/(GLfloat)h, 0.0, 1.0);*/
}

void keyboard(unsigned char key, int x, int y)
{
   switch (key) {
      case 'd':
	  case 'D':
		  /* do something here*/
		 
		  glutPostRedisplay();
		  break;
	  case 'a':
	  case 'A':
		  /* do something here*/
		  // your commands
		
		  glutPostRedisplay();
         break;
      case 27:  /*  Escape key  */
         exit(0);
         break;
      default:
         break;
   }
}

void mouseClick(int btn, int state, int x, int y)
{
	if (state == GLUT_UP)
	{
		x_before = 0;
		y_before = 0;
	}
	if (btn == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
	{
		bubbleGame.findBubble(x, CANVASH- y);
		//cout << "hhhh" << endl;
		glutPostRedisplay();
		//display();
	}
	bubbleGame.clearSlected();
}

void mouseMove(int x, int y) {
	int y_convert = CANVASH - y;
	int x_move = 0;
	int y_move = 0;

	if (x_before == 0 && y_before == 0)
	{
		x_move = 0;
		y_move = 0;
	}
	else
	{
		x_move = x - x_before;
		y_move = y_convert - y_before;
	}
	x_before = x;
	y_before = y_convert;
	bubbleGame.findCluster(x, y_convert, x_move, y_move);
	//bubbleGame.moveBubble(x, y_convert);
	glutPostRedisplay();
	//display();
	//bubbleGame.clustering_spliting();
	//bubbleGame.drawCurve();
	
	//bubbleGame.movedDrawCurve(x, y_convert);
}



void timer(int value)
{
	// do something here if you have enabled a timer

	glutPostRedisplay();
	glutTimerFunc(700, timer, 0);
}

void drawPoint()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glPointSize(5.0f);
	glBegin(GL_POINTS);
	glVertex2f(0.0f, 0.0f);
	glVertex2f(0.5f, 0.5f);
	glEnd();
	glFlush();
}

/*  Main Loop
 *  Open window with initial window size, title bar, 
 *  RGBA display mode, and handle input events.
 */
int main(int argc, char** argv)
{
   glutInit(&argc, argv);
   glutInitDisplayMode(GLUT_DOUBLE|GLUT_RGB|GLUT_DEPTH);
   glutInitWindowSize (800, 600);
   glutCreateWindow (argv[0]);
   //drawPoint();
   init();
   glutReshapeFunc (reshape);
   glutKeyboardFunc (keyboard);
   glutMouseFunc (mouseClick);
   glutDisplayFunc(display);
   glutMotionFunc(mouseMove);
//   glutTimerFunc (700, timer, 0);
   glutMainLoop();
   return 0;
}
