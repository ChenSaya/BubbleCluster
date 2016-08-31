#include "Bubble.h"
#include <GL/glut.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <cmath>

/* hints
// draw fonts
void RenderString(float x, float y, void *font, const char* string, RGB const& rgb)
{  
  char *c;

  glColor3f(rgb.r, rgb.g, rgb.b); 
  glRasterPos2f(x, y);

  glutBitmapString(font, string);
}
// pause the game
//Set handler functions for drawing, keypresses, and window resizes                                                                                              
if(!paused)
{
	glutDisplayFunc(drawScene);       //drawScene draws everything does the real work
	glutTimerFunc(10, update, 0); //Add a timer                                                                                                                     
}
*/






Icon::~Icon()
{
	
}

void Icon::draw()
{
	float x1 = x_, x2 = x_ + w_;
	float y1 = y_, y2 = y_ + h_;
	// draw line
	glColor4f(0, 0, 0, 1);
	glLineWidth(2);
	glBegin (GL_LINES);
	glVertex2f(x1,y1);
	glVertex2f(x1,y2);
	glVertex2f(x1,y1);
	glVertex2f(x2,y1);
	glVertex2f(x1,y2);
	glVertex2f(x2,y2);
	glVertex2f(x2,y1);
	glVertex2f(x2,y2);
	glEnd();
	// draw quad
	glColor4f(1, 1, 0, 0.75);
	glBegin (GL_QUADS);
	glVertex2f(x1,y1);
	glVertex2f(x2,y1);
	glVertex2f(x2,y2);
	glVertex2f(x1,y2);
	glEnd();
}



void Icon::moveIcon(int x, int y, int x_move, int y_move)
{
	float x1 = x_, x2 = x_ + w_;
	float y1 = y_, y2 = y_ + h_;
	if (x == 0 && y == 0)
	{
		iconSlected = false;
	}

	if (x1 <x && x<x2 && y>y1 && y<y2 && iconSlected == false)
	{
		iconSlected = true;
	}
	if(iconSlected)
	{
		x_ = x_ + x_move;
		y_ = y_ + y_move;
	}	
}







BubbleCluster::BubbleCluster(vector<Icon*> icons)
{
	this->icons_ = icons;
}

BubbleCluster::BubbleCluster(Icon* icon)
{
	this->addIcon(icon);
}



BubbleCluster::BubbleCluster()
{
	
}

BubbleCluster::~BubbleCluster()
{
	// release icons
	vector<Icon*>::iterator icon_iter;
	for (icon_iter = icons_.begin(); icon_iter != icons_.end(); ++icon_iter)
	{
		delete *icon_iter;
	}
	this->icons_.clear(); // clear content
}


void BubbleCluster::draw()
{
	// draw cluster boundary here

	
	//drawCurve();
	/*int count;

	float x = 0, y=0;
	int sections = 200;
	GLfloat TWOPI = 2.0f * 3.14159f;
	

	
	for (icon_iter = icons_.begin(); icon_iter != icons_.end(); ++icon_iter)
	{
		(*icon_iter)->return_position(&x, &y);
		glBegin(GL_LINE_STRIP);
		for (count = 0; count <= sections; count++)
		{
			glVertex2f(x+12 + Radius*cos(count*TWOPI / sections), y+12 + Radius*sin(count*TWOPI / sections));
		}
		glEnd();

	}
	*/
	vector<Icon*>::iterator icon_iter;
	// draw icons
	for (icon_iter = icons_.begin(); icon_iter != icons_.end(); ++icon_iter)
	{
		(*icon_iter)->draw();
	}
}

/*void BubbleCluster::return_icon(vector<Icon*> icons)
{
	icons_ = icons;
}*/

void BubbleCluster::findIcon(int x, int y, int x_move, int y_move)
{
	vector<Icon*>::iterator icon_iter;
	for (icon_iter = icons_.begin(); icon_iter != icons_.end(); ++icon_iter)
	{
		(*icon_iter)->moveIcon(x, y, x_move, y_move);
	}
}
void BubbleCluster::clearSlected()
{
	vector<Icon*>::iterator icon_iter;
	for (icon_iter = icons_.begin(); icon_iter != icons_.end(); ++icon_iter)
	{
		(*icon_iter)->clearSlected();
	}
}

void BubbleCluster::clusterWeight(int x, int y)
{
	int cluster_size = icons_.size();
	int x_icon, y_icon;
	float h0 = hmax;
	float* point_to_icon_distance = new float[cluster_size];
	for (int i = 0; i < cluster_size; i++)
	{
		icons_[i]->getIconPosition(&x_icon, &y_icon);
		point_to_icon_distance[i] = radialQuadraticFunction(x, x_icon, y, y_icon);
		h0 = h0 + point_to_icon_distance[i];
	}
	float* point_to_icon_weight = new float[cluster_size];
	float fx = 0;
	for (int j = 0; j < cluster_size; j++)
	{
		point_to_icon_weight[j] = hmax / h0;
		fx = fx + point_to_icon_weight[j] * point_to_icon_distance[j];
	}
	delete[] point_to_icon_distance;
	delete[] point_to_icon_weight;
	setClusterFx(fx);
}
void BubbleCluster::boundingBox()
{
	int boxxmin, boxxmax, boxymin, boxymax;
	int xmin, ymin, xmax, ymax, x, y;
	icons_[0]->getIconPosition(&xmin, &ymin);
	icons_[0]->getIconPosition(&xmax, &ymax);
	vector<Icon*>::iterator icon_iter;
	for (icon_iter = icons_.begin(); icon_iter != icons_.end(); ++icon_iter)
	{
		(*icon_iter)->getIconPosition(&x, &y);
		if (x >= xmax)
		{
			xmax = x;
		}
		if (y >= ymax)
		{
			ymax = y;
		}
		if (x <= xmin)
		{
			xmin = x;
		}
		if (y <= ymin)
		{
			ymin = y;
		}
	}
	/*if (xmin - CLUSTERBOUND < 0) boxxmin = 0;
	if (ymin - CLUSTERBOUND < 0) boxymin = 0;
	if (xmax + CLUSTERBOUND > CANVASW) boxxmax = CANVASW;
	if (ymax + CLUSTERBOUND > CANVASH) boxymax = CANVASH;
	else
	{*/
	boxxmin = xmin - CLUSTERBOUND;// +ICONSIZE / 2;
	boxymin = ymin - CLUSTERBOUND;// +ICONSIZE / 2;
	boxxmax = xmax + CLUSTERBOUND;// +ICONSIZE / 2;
	boxymax = ymax + CLUSTERBOUND;// +ICONSIZE / 2;
	
	setBoundingBox(boxxmin, boxymin, boxxmax, boxymax);
}
bool BubbleCluster::clusterSlected(int x, int y)
{
	int iconxmin, iconymin, iconxmax, iconymax;
	vector<Icon*>::iterator icon_iter;
	for (icon_iter = icons_.begin(); icon_iter != icons_.end(); ++icon_iter)
	{
		(*icon_iter)->getIconPosition(&iconxmin, &iconymin);
		iconxmax = iconxmin + ICONSIZE;
		iconymax = iconymin + ICONSIZE;
		if (x >= iconxmin&&x <= iconxmax&&y >= iconymin&&y <= iconymax)
		{
			return true;
		}
	}
	return false;
}
void BubbleCluster::pointToClusterValue(int x, int y)
{
	int cluster_size = icons_.size();
	int x_icon, y_icon;
	float h0 = hmax;
	float* point_to_icon_distance = new float[cluster_size];
	for (int i = 0; i < cluster_size; i++)
	{
		icons_[i]->getIconPosition(&x_icon, &y_icon);
		point_to_icon_distance[i] = radialQuadraticFunction(x, x_icon, y, y_icon);
		h0 = h0 + point_to_icon_distance[i];
	}
	float* point_to_icon_weight = new float[cluster_size];
	float fx = 0;
	for (int j = 0; j < cluster_size; j++)
	{
		point_to_icon_weight[j] = hmax / h0;
		fx = fx + point_to_icon_weight[j] * point_to_icon_distance[j];
	}
	delete[] point_to_icon_distance;
	delete[] point_to_icon_weight;
	setClusterFx(fx);
}

void BubbleCluster::triangulateInit(triangulateio * io)
{
	io->pointlist = 0;
	io->pointattributelist = 0;
	io->pointmarkerlist = 0;
	io->numberofpoints = 0;
	io->numberofpointattributes = 0;

	io->trianglelist = 0;
	io->triangleattributelist = 0;
	io->trianglearealist = 0;
	io->neighborlist = 0;
	io->numberoftriangles = 0;
	io->numberofcorners = 0;
	io->numberoftriangleattributes = 0;

	io->segmentlist = 0;
	io->segmentmarkerlist = 0;
	io->numberofsegments = 0;

	io->holelist = 0;
	io->numberofholes = 0;

	io->regionlist = 0;
	io->numberofregions = 0;

	io->edgelist = 0;
	io->edgemarkerlist = 0;
	io->normlist = 0;
	io->numberofedges = 0;
}

void BubbleCluster::triangulateT()
{
	struct triangulateio in, mid, vorout;
	printf("Initial triangulation:\n\n");
	triangulateInit(&in);
	triangulateInit(&mid);
	triangulateInit(&vorout);
	/* Define input points. */


	in.numberofpoints = icons_.size();
	cout <<"bubble size is "<< icons_.size() << endl;
	in.numberofsegments = icons_.size();
	in.pointlist = (REAL *)malloc(in.numberofpoints * 2 * sizeof(REAL));
	in.pointmarkerlist = (int *)malloc(in.numberofpoints * sizeof(int));
	in.segmentlist = (int*)malloc(in.numberofpoints * 2 * sizeof(int));
	in.segmentmarkerlist = (int *)malloc(in.numberofpoints * sizeof(int));
	// input points
	for (int i = 0; i < in.numberofpoints; i++)
	{
		int x, y;
		icons_[i]->getIconPosition(&x, &y);
		cout << "x is " << x << " | y is " << y << endl;
		in.pointlist[i << 1] = x;
		in.pointlist[(i << 1) + 1] = y;
		in.pointmarkerlist[i] = 0;
	//	in.segmentlist[i << 1] = i;
	//	in.segmentlist[(i << 1) + 1] = (i + 1) % (in.numberofpoints);
	//	in.segmentmarkerlist[i] = 1;
	}
	triangulate("pczAen", &in, &mid, &vorout);
	out_tri_ = mid;
	trianguleReport(&mid, 1, 1, 1, 1, 1, 0);
	cout << "end trianglate " << endl;

	//triangulate("prazBP", &mid, &out, (struct triangulateio *) NULL);
	//addMesh(out);

	if (vorout.pointlist != 0)			free(vorout.pointlist);
	if (vorout.pointattributelist != 0)	free(vorout.pointattributelist);
	if (vorout.edgelist != 0)			free(vorout.edgelist);
	if (vorout.normlist != 0)			free(vorout.normlist);
}

void BubbleCluster::trianguleReport(triangulateio * io, int markers, int reporttriangles, int reportneighbors, int reportsegments, int reportedges, int reportnorms)
{
	int i, j;

	for (i = 0; i < io->numberofpoints; i++) {
		printf("Point %4d:", i);
		for (j = 0; j < 2; j++) {
			printf("  %.6g", io->pointlist[i * 2 + j]);
		}
		if (io->numberofpointattributes > 0) {
			printf("   attributes");
		}
		for (j = 0; j < io->numberofpointattributes; j++) {
			printf("  %.6g",
				io->pointattributelist[i * io->numberofpointattributes + j]);
		}
		if (markers) {
			printf("   marker %d\n", io->pointmarkerlist[i]);
		}
		else {
			printf("\n");
		}
	}
	printf("\n");

	if (reporttriangles || reportneighbors) {
		for (i = 0; i < io->numberoftriangles; i++) {
			if (reporttriangles) {
				printf("Triangle %4d points:", i);
				for (j = 0; j < io->numberofcorners; j++) {
					printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
				}
				if (io->numberoftriangleattributes > 0) {
					printf("   attributes");
				}
				for (j = 0; j < io->numberoftriangleattributes; j++) {
					printf("  %.6g", io->triangleattributelist[i *
						io->numberoftriangleattributes + j]);
				}
				printf("\n");
			}
			if (reportneighbors) {
				printf("Triangle %4d neighbors:", i);
				for (j = 0; j < 3; j++) {
					printf("  %4d", io->neighborlist[i * 3 + j]);
				}
				printf("\n");
			}
		}
		printf("\n");
	}

	if (reportsegments) {
		for (i = 0; i < io->numberofsegments; i++) {
			printf("Segment %4d points:", i);
			for (j = 0; j < 2; j++) {
				printf("  %4d", io->segmentlist[i * 2 + j]);
			}
			if (markers) {
				printf("   marker %d\n", io->segmentmarkerlist[i]);
			}
			else {
				printf("\n");
			}
		}
		printf("\n");
	}

	if (reportedges) {
		for (i = 0; i < io->numberofedges; i++) {
			printf("Edge %4d points:", i);
			for (j = 0; j < 2; j++) {
				printf("  %4d", io->edgelist[i * 2 + j]);
			}
			if (reportnorms && (io->edgelist[i * 2 + 1] == -1)) {
				for (j = 0; j < 2; j++) {
					printf("  %.6g", io->normlist[i * 2 + j]);
				}
			}
			if (markers) {
				printf("   marker %d\n", io->edgemarkerlist[i]);
			}
			else {
				printf("\n");
			}
		}
		printf("\n");
	}
}




void BubbleCluster::spreadingTwo()
{
	int icon_x1, icon_y1, icon_x2, icon_y2;
	icons_[0]->getIconPosition(&icon_x1, &icon_y1);
	icons_[1]->getIconPosition(&icon_x2, &icon_y2);
	if (icon_x1 <= icon_x2 && icon_y1 < icon_y2)
	{
		icon_x1 = icon_x1 - SPREAD_DISTANCE_TWO;
		icon_x2 = icon_x2 + SPREAD_DISTANCE_TWO;

		icon_y1 = icon_y1 - SPREAD_DISTANCE_TWO;
		icon_y2 = icon_y2 + SPREAD_DISTANCE_TWO;
	}
	if (icon_x1 <= icon_x2 && icon_y1 > icon_y2)
	{
		icon_x1 = icon_x1 - SPREAD_DISTANCE_TWO;
		icon_x2 = icon_x2 + SPREAD_DISTANCE_TWO;

		icon_y1 = icon_y1 + SPREAD_DISTANCE_TWO;
		icon_y2 = icon_y2 - SPREAD_DISTANCE_TWO;
	}
	if (icon_x1 > icon_x2 && icon_y1 < icon_y2)
	{
		icon_x1 = icon_x1 + SPREAD_DISTANCE_TWO;
		icon_x2 = icon_x2 - SPREAD_DISTANCE_TWO;

		icon_y1 = icon_y1 - SPREAD_DISTANCE_TWO;
		icon_y2 = icon_y2 + SPREAD_DISTANCE_TWO;
	}
	if (icon_x1 > icon_x2 && icon_y1 > icon_y2)
	{
		icon_x1 = icon_x1 + SPREAD_DISTANCE_TWO;
		icon_x2 = icon_x2 - SPREAD_DISTANCE_TWO;

		icon_y1 = icon_y1 + SPREAD_DISTANCE_TWO;
		icon_y2 = icon_y2 - SPREAD_DISTANCE_TWO;
	}
	else
	{
		icon_x1 = icon_x1 + ICONSIZE;
		icon_x2 = icon_x2 - ICONSIZE;

		icon_y1 = icon_y1 + ICONSIZE;
		icon_y2 = icon_y2 - ICONSIZE;
	}
	icons_[0]->setIconPosition(icon_x1, icon_y1);
	icons_[1]->setIconPosition(icon_x2, icon_y2);
}

void BubbleCluster::spreading()
{
	int m = out_tri_.numberofedges;
	int n = out_tri_.numberofpoints;

	for (int i = 0; i < out_tri_.numberofpoints; i++)
	{
		//cout << "point x is: " << out_tri_.pointlist[2 * i];
		//cout << "   |    point y is: " << out_tri_.pointlist[2 * i + 1] << endl;
	}

	bool overlapped = false;
	for (int i = 0; i < n; i++)
	{
		overlapped = false;
		int point_x_i = out_tri_.pointlist[2 * i];
		int point_y_i = out_tri_.pointlist[2 * i + 1];
		for (int j = i + 1; j < n && overlapped == false; j++)
		{
			int point_x_j = out_tri_.pointlist[2 * j];
			int point_y_j = out_tri_.pointlist[2 * j + 1];
			//cout << "befor relative : x is : " << (point_x_i - point_x_j) << " | y is " << (point_y_i - point_y_j) << endl;
			if (abs(point_x_i - point_x_j) < ICONSIZE && abs(point_y_i - point_y_j) < ICONSIZE)
			{
				overlapped = true;
				//cout << "overlap" << endl;
				if (point_x_i < point_x_j && point_y_i < point_y_j)
				{
					point_x_i = point_x_i - SPREAD_DISTANCE;
					point_x_j = point_x_j + SPREAD_DISTANCE;

					point_y_i = point_y_i - SPREAD_DISTANCE;
					point_y_j = point_y_j + SPREAD_DISTANCE;
				}
				if (point_x_i < point_x_j && point_y_i > point_y_j)
				{
					point_x_i = point_x_i - SPREAD_DISTANCE;
					point_x_j = point_x_j + SPREAD_DISTANCE;

					point_y_i = point_y_i + SPREAD_DISTANCE;
					point_y_j = point_y_j - SPREAD_DISTANCE;
				}
				if (point_x_i > point_x_j && point_y_i < point_y_j)
				{
					point_x_i = point_x_i + SPREAD_DISTANCE;
					point_x_j = point_x_j - SPREAD_DISTANCE;

					point_y_i = point_y_i - SPREAD_DISTANCE;
					point_y_j = point_y_j + SPREAD_DISTANCE;
				}

				if (point_x_i > point_x_j && point_y_i > point_y_j)
				{
					point_x_i = point_x_i + SPREAD_DISTANCE;
					point_x_j = point_x_j - SPREAD_DISTANCE;

					point_y_i = point_y_i + SPREAD_DISTANCE;
					point_y_j = point_y_j - SPREAD_DISTANCE;
				}
				else
				{
					point_x_i = point_x_i + ICONSIZE;
					point_x_j = point_x_j - ICONSIZE;

					point_y_i = point_y_i + ICONSIZE;
					point_y_j = point_y_j - ICONSIZE;
				}

				out_tri_.pointlist[i * 2] = point_x_i;
				out_tri_.pointlist[i * 2 + 1] = point_y_i;

				out_tri_.pointlist[j * 2] = point_x_j;
				out_tri_.pointlist[j * 2 + 1] = point_y_j;
			}
		}
		if (overlapped)
		{
			//cout << "i is: " << i << endl;
			i = 0;
			//cout << "8888888" << endl;
		}
	}
	
	/*for (i = 0; i < clusters.size(); i++)
	{
		icons_1 = clusters[i]->return_icon();
		ifMerged = false;
		for (j = i + 1; j < clusters.size() && ifMerged == false; j++)
		{
			icons_2 = clusters[j]->return_icon();
			if (ifClusterMerge(icons_1, icons_2))
			{
				clusters[i]->addIcon(icons_2[0]);
				removeCluster(&clusters, clusters[j]);
				//cluster_size = clusters.size();
				ifMerged = true;
				cout << "cluster size is " << clusters.size() << endl;
			}
		}
		if (ifMerged)
		{
			i = 0;//if merge

		}
	}*/
	/*
	bool overlapped = false;
	for (int i = 0; i < out_tri_.numberofpoints; i++)
	{
		overlapped = false;
		int point_index_left = out_tri_.edgelist[i *2];
		int point_index_right = out_tri_.edgelist[i * 2 + 1];
		cout << "left point index is: " << point_index_left;
		cout << "right point index is: " << point_index_right << endl;
		int left_x = out_tri_.pointlist[point_index_left * 2];
		int left_y = out_tri_.pointlist[point_index_left * 2 + 1];

		int right_x = out_tri_.pointlist[point_index_right * 2];
		int right_y = out_tri_.pointlist[point_index_right * 2 + 1];
		cout << "befor relative : x is : " << (left_x - right_x) << " | y is " << (left_y - right_y) << endl;
		if (abs(left_x - right_x) < ICONSIZE && abs(left_y - right_y) < ICONSIZE && overlapped == false)
		{
			overlapped = true;
			cout << "overlap" << endl;
			left_x = left_x + ICONSIZE * 2;
			left_y = left_y + ICONSIZE * 2;
			out_tri_.pointlist[point_index_left * 2] = left_x;
			out_tri_.pointlist[point_index_left * 2 + 1] = left_y;
		}
		if (overlapped == true)
		{
			i = 0;
		}
	}
	*/



	//test
	/*
	for (int i = 0; i < n; i++)
	{
		cout << "after point x is: " << out_tri_.pointlist[2 * i];
		cout << "   |   after point y is: " << out_tri_.pointlist[2 * i + 1] << endl;
	}
	*/


	//cout << "n is: " << n << endl;
	//cout << "m is: " << m << endl;
	MatrixXf b_x(m + 1, 1);
	MatrixXf b_y(m + 1, 1);

	for (int i = 0; i < m; i++)
	{
		int point_index_left = out_tri_.edgelist[i * 2];
		int point_index_right = out_tri_.edgelist[i * 2 + 1];
		int relative_position_x = out_tri_.pointlist[point_index_left * 2] - out_tri_.pointlist[point_index_right * 2];
		int relative_position_y = out_tri_.pointlist[point_index_left * 2 + 1] - out_tri_.pointlist[point_index_right * 2 + 1];
		//cout << "relative x is: " << relative_position_x << "   ralative y is: " << relative_position_y << endl;
		b_x(i, 0) = relative_position_x;
		b_y(i, 0) = relative_position_y;
		
	}

	int center_x;
	int center_y;
	for (int i = 0; i < n; i++)
	{
		center_x += (out_tri_.pointlist[i * 2])/n;
		center_y += (out_tri_.pointlist[i * 2 + 1]) / n;
	}
	b_x(m, 0) = center_x;
	b_y(m, 0) = center_y;
	//cout << "center x is : " << center_x << "  |  center y is: " << center_y << endl;



	MatrixXf A(m + 1, n);

	for (int i = 0; i < m + 1; i++)
	{
		for (int j = 0; j < n; j++)
		{
			A(i, j) = 0;
		}
	}
	
	for (int i = 0; i < m; i++)
	{
		A(i, out_tri_.edgelist[i * 2]) = 1;
		A(i, out_tri_.edgelist[i * 2 + 1]) = -1;
	}
	for (int i = 0; i < n; i++)
	{
		A(m, i) = 1.0 / n;
		//cout << "the last row is: "<< A(m, i) << endl;
	}
	for (int i = 0; i < m + 1; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << A(i, j) << " ";
		}
		cout << endl;
	}


	// calculate v'

	MatrixXf v_x(n, 1);
	MatrixXf v_y(n, 1);


	MatrixXf A_T = (A.transpose() * A).inverse() * A.transpose();
	//cout << "A_T" << endl<<A_T << endl;
	//cout << "b_x " << endl << b_x << endl;
	//cout << "b_y " << endl << b_y << endl;
	v_x = A_T * b_x;
	v_y = A_T * b_y;


	//cout << "v_x is: " << endl << v_x << endl;
	//cout << "v_y is: " << endl << v_y << endl;



	for (int i = 0; i < icons_.size(); i++)
	{
		int x, y;
		icons_[i]->setIconPosition(v_x(i, 0), v_y(i, 0));
	}
	
}









BubbleGameEngine::BubbleGameEngine()
{
	vector<BubbleCluster*>::iterator cluster_iter;
	for (cluster_iter = bubble_clusters_.begin(); cluster_iter != bubble_clusters_.end(); ++cluster_iter)
	{
		*cluster_iter;
	}
}

BubbleGameEngine::~BubbleGameEngine()
{
	// release bubble clusters
	vector<BubbleCluster*>::iterator cluster_iter;
	for (cluster_iter = bubble_clusters_.begin(); cluster_iter != bubble_clusters_.end(); ++cluster_iter)
	{
		delete *cluster_iter;
	}
}

void BubbleGameEngine::init()
{
	srand(time(NULL));

	//int K = 10; 
	for (int i = 0; i < K; ++i)
	{
		float x = rand() % CANVASW;
		float y = rand() % CANVASH;
		if (x + ICONSIZE > CANVASW) x = CANVASW - ICONSIZE - 5;
		if (y + ICONSIZE > CANVASH) y = CANVASH - ICONSIZE - 5;
		Icon* icon = new Icon(x, y, ICONSIZE, ICONSIZE);
		icon->setGroupLable(i+1);
		//cout << " "<<icon->returnGroupLable()<<" " << endl;
		BubbleCluster* bb_cluster = new BubbleCluster(icon);
		this->bubble_clusters_.push_back(bb_cluster);
		bb_cluster->setGroupLable(i + 1);
	}
	clusteringInitial(bubble_clusters_);
	//drawBoundingBox();
}

void BubbleGameEngine::draw()
{
	
	vector<BubbleCluster*>::iterator cluster_iter;
	for (cluster_iter = bubble_clusters_.begin(); cluster_iter != bubble_clusters_.end(); ++cluster_iter)
	{
		(*cluster_iter)->draw();
	}
	clustering_spliting();
	//drawCurve();
	drawCurve2();
}

void BubbleGameEngine::findCluster(int x, int y, int x_move, int y_move)
{
	vector<BubbleCluster*>::iterator cluster_iter;
	for (cluster_iter = bubble_clusters_.begin(); cluster_iter != bubble_clusters_.end(); ++cluster_iter)
	{
		(*cluster_iter)->findIcon(x, y, x_move, y_move);
	}
}

void BubbleGameEngine::clearSlected()
{
	vector<BubbleCluster*>::iterator cluster_iter;
	for (cluster_iter = bubble_clusters_.begin(); cluster_iter != bubble_clusters_.end(); ++cluster_iter)
	{
		(*cluster_iter)->clearSlected();
	}
}

bool BubbleGameEngine::ifClusterMerge(vector<Icon*> icons_i, vector<Icon*> icons_j)
{
	int m, n;
	for (m = 0; m < icons_i.size(); m++)
	{
		for (n = 0; n < icons_j.size(); n++)
		{
			if (calculateDistance(icons_i[m], icons_j[n]) < SMALL_THREASHOLD)
			{
				return true;
			}
		}
	}
	return false;

}

// need to complete
vector<BubbleCluster*> BubbleGameEngine::clusteringInitial(vector<BubbleCluster*> clusters)
{
	//vector<BubbleCluster*> tempclusters;
	int i, j;
	bool ifMerged = false;
	vector<Icon*> icons_1, icons_2;
	//int cluster_size = clusters.size();
	for (i = 0; i < clusters.size(); i++)
	{
		icons_1 = clusters[i]->return_icon();
		ifMerged = false;
		for (j = i + 1; j < clusters.size() && ifMerged == false; j++)
		{
			icons_2 = clusters[j]->return_icon();
			if (ifClusterMerge(icons_1, icons_2))
			{
				clusters[i]->addIcon(icons_2[0]);
				removeCluster(&clusters, clusters[j]);
				//cluster_size = clusters.size();
				ifMerged = true;
				//cout << "cluster size is " << clusters.size() << endl;
			}
		}
		if (ifMerged)
		{
			i = 0;//if merge

		}
	}
	return clusters;
}

void BubbleGameEngine::clustering_spliting()
{

	/*	vector<BubbleCluster*>::iterator cluster_iter1;
		vector<BubbleCluster*>::iterator cluster_iter2;
		vector<Icon*>::iterator icon_iter1;
		vector<Icon*>::iterator icon_iter2;
		
		
		for (cluster_iter1 = bubble_clusters_.begin(); cluster_iter1 != bubble_clusters_.end(); ++cluster_iter1)
		{
			vector<Icon*> icons_1 = (*cluster_iter1)->return_icon();	
			for (icon_iter1 = icons_1.begin(); icon_iter1 != icons_1.end(); ++icon_iter1)
			{
				for (cluster_iter2 = bubble_clusters_.begin(); cluster_iter2 != bubble_clusters_.end(); ++cluster_iter2)
				{
					vector<Icon*> icons_2 = (*cluster_iter2)->return_icon();
					for (icon_iter2 = icons_2.begin(); icon_iter2 != icons_2.end(); ++icon_iter2)
					{
						if (calculateDistance(*icon_iter1, *icon_iter2) < threashold(*icon_iter1, *icon_iter2))
						{
							(*icon_iter2)->setGroupLable((*icon_iter2)->returnGroupLable());
							(*cluster_iter1)->addIcon(*icon_iter2);
							(*cluster_iter2)->removeIcon(*icon_iter2);
							int cluster_length = (*cluster_iter2)->return_icon().size();
							cout << "current bubble size" << bubble_clusters_.size() << endl;
							int length = bubble_clusters_.size();
							if (cluster_length == 0)
							{
								removeCluster(*cluster_iter2);
								cout << "after remove cluster" << bubble_clusters_.size() << endl;
								cout << "Cluster " << (*cluster_iter2)->returnGroupLable() << " has been removed" << endl;
								return;
							}
						}
					}
				}
			}
			
		}*/

	int cluster_size = bubble_clusters_.size();
	//cout << "cluster number:" << cluster_size << endl;

	for (int i = 0; i < cluster_size; i++)
	{
		for (int j = 0; j < cluster_size; j++)
		{			
			if (i == j) // this means we have selected the same cluster, then spliting
			{
				vector<Icon*> cluster_split_i = (bubble_clusters_[i])->return_icon();
				vector<Icon*> cluster_split_j = (bubble_clusters_[j])->return_icon();
				int len_split_i = cluster_split_i.size(); // the ith cluster length
				int len_split_j = cluster_split_j.size(); // the jth cluster length
				for (int k_split = 0; k_split < len_split_i; k_split++)
				{
					for (int l_split = 0; l_split < len_split_j; l_split++)
					{
						if (calculateDistance(cluster_split_i[k_split], cluster_split_j[l_split]) > LARGE_THREASHOLD)
						{
							BubbleCluster* split_new_cluster = new BubbleCluster(cluster_split_j[l_split]);
							bubble_clusters_.push_back(split_new_cluster);
							bubble_clusters_[j]->removeIcon(cluster_split_j[l_split]);
							cluster_size = bubble_clusters_.size();
							//cout << "Cluster size" << cluster_size << endl;
							return;
						}
					}
				}
			}
			else {
				//j = i + 1;
				int group_i_label;
				vector<Icon*> cluster_i = (bubble_clusters_[i])->return_icon();
				vector<Icon*> cluster_j = (bubble_clusters_[j])->return_icon();
				int len_i = cluster_i.size(); // the ith cluster length
				int len_j = cluster_j.size(); // the jth cluster length
				cluster_i[0]->getGroupLable(&group_i_label);
				for (int k = 0; k < len_i; k++)
				{
					for (int l = 0; l < len_j; l++)
					{
						if (calculateDistance(cluster_i[k], cluster_j[l]) < SMALL_THREASHOLD)
						{
							cluster_j[l]->setGroupLable(group_i_label); // set the jth group label equal to ith group
							bubble_clusters_[i]->addIcon(cluster_j[l]); // move jth cluster's icons into ith cluster
							bubble_clusters_[j]->removeIcon(cluster_j[l]);
							int leni = bubble_clusters_[j]->return_icon().size();
							int leng = bubble_clusters_.size();
							if (bubble_clusters_[j]->return_icon().size() == 0)
							{
								removeCluster(bubble_clusters_[j]); // remove jth cluster
								cluster_size = bubble_clusters_.size();
								cout << "Cluster size" << cluster_size << endl;
								//cout << "Cluster " << bubble_clusters_[j]->returnGroupLable() << " has been removed" << endl;
								return;
							}


						}
					}

				}
			}
		}
	}	
}

int BubbleGameEngine::calculateDistance(Icon* icon_i, Icon* icon_j)
{
	// x,y  means the icon center coordinate, x_ and y_ means other coordinate: mouse or something
	int x_i, y_i, x_j, y_j;
	icon_i->getIconPosition(&x_i, &y_i);
	icon_j->getIconPosition(&x_j, &y_j);
	int distance = sqrt(pow(x_i - x_j, 2) + pow(y_i - y_j, 2));
	return distance;
}


void BubbleGameEngine::drawCurve()
{
	//int count = 0;
	for (int x = 0; x < CANVASW; x++)
	{
		for (int y = 0; y < CANVASH; y++)
		{
			//int bubble_size = bubble_clusters_.size();
			float sum_cluster_fx = 0.0;
			float fx_i, fx_max;
			vector<BubbleCluster*>::iterator bubble_iter;
			for (bubble_iter = bubble_clusters_.begin(); bubble_iter != bubble_clusters_.end(); ++bubble_iter)
			{
				(*bubble_iter)->clusterWeight(x, y);
				(*bubble_iter)->getClusterFx(&fx_i);
				sum_cluster_fx = sum_cluster_fx + fx_i;
			}
			for (bubble_iter = bubble_clusters_.begin(); bubble_iter != bubble_clusters_.end(); ++bubble_iter)
			{
				(*bubble_iter)->getClusterFx(&fx_max);
				
				if (fx_max*(1 + NEGATIVEAFFECT) - NEGATIVEAFFECT*sum_cluster_fx >= 1)
				{
					glPointSize(2);
					glColor3f(0, 1, 1);
					glBegin(GL_POINTS);
					
					if ((1 < fx_max*(1 + NEGATIVEAFFECT) - NEGATIVEAFFECT*sum_cluster_fx) && (fx_max*(1 + NEGATIVEAFFECT) - NEGATIVEAFFECT*sum_cluster_fx < 1.07))
					{
						float x_ = x + ICONSIZE/2;
						float y_ = y + ICONSIZE/2;
						glVertex2f(x_, y_);
						//count++;
					}
					glEnd();	
					//return;
				}
			}
			/*for (int i = 0; i < bubble_size; i++)
			{
				bubble_clusters_[i]->clusterWeight(x, y);
				float fx_i = bubble_clusters_[i]->returnClusterFx();
				sum_cluster_fx = sum_cluster_fx + fx_i;
			}
			for (int i = 0; i < bubble_size; i++)
			{
				float fx_max = bubble_clusters_[i]->returnClusterFx();
				if (2 * fx_max - sum_cluster_fx >= 1)
				{
					glPointSize(3);
					glBegin(GL_POINTS);
					if (1 < fx_max && fx_max < 1.03)
					{
						glVertex2f(x, y);
					}
					glEnd();
				}
			}*/

		}
	}
}

void BubbleGameEngine::drawCurve2()
{
	int xmin, ymin, xmax, ymax;
	vector<BubbleCluster*>::iterator cluster_iter;
	for (cluster_iter = bubble_clusters_.begin(); cluster_iter != bubble_clusters_.end(); ++cluster_iter)
	{
		(*cluster_iter)->boundingBox();
		(*cluster_iter)->getBoundingBox(&xmin, &ymin, &xmax, &ymax);
		for (int x = xmin; x < xmax; x++)
		{
			for (int y = ymin; y < ymax; y++)
			{
				//int bubble_size = bubble_clusters_.size();
				float sum_cluster_fx = 0.0;
				float fx_i, fx_max;
				vector<BubbleCluster*>::iterator bubble_iter;
				for (bubble_iter = bubble_clusters_.begin(); bubble_iter != bubble_clusters_.end(); ++bubble_iter)
				{
					(*bubble_iter)->clusterWeight(x, y);
					(*bubble_iter)->getClusterFx(&fx_i);
					sum_cluster_fx = sum_cluster_fx + fx_i;
				}
				for (bubble_iter = bubble_clusters_.begin(); bubble_iter != bubble_clusters_.end(); ++bubble_iter)
				{
					(*bubble_iter)->getClusterFx(&fx_max);

					if (fx_max*(1 + NEGATIVEAFFECT) - NEGATIVEAFFECT*sum_cluster_fx >= 1)
					{
						glPointSize(2);
						glColor3f(0, 1, 1);
						glBegin(GL_POINTS);

						if ((1 < fx_max*(1 + NEGATIVEAFFECT) - NEGATIVEAFFECT*sum_cluster_fx) && (fx_max*(1 + NEGATIVEAFFECT) - NEGATIVEAFFECT*sum_cluster_fx < 1.03))
						{
							float x_ = x + ICONSIZE / 2;
							float y_ = y + ICONSIZE / 2;
							glVertex2f(x_, y_);
							//count++;
						}
						glEnd();
						//return;
					}
				}
			}
		}

	}
}



void BubbleGameEngine::drawBoundingBox()
{
	int x_min, y_min, x_max, y_max;
	vector<BubbleCluster*>::iterator cluster_iter;
	for (cluster_iter = bubble_clusters_.begin(); cluster_iter != bubble_clusters_.end(); ++cluster_iter)
	{
		(*cluster_iter)->boundingBox();
		(*cluster_iter)->getBoundingBox(&x_min, &y_min, &x_max, &y_max);
		for (int x = 0; x < CANVASW; x++)
		{
			for (int y = 0; y < CANVASH; y++)
			{

				glPointSize(5);
				glColor4f(0, 0, 0, 1);
				glBegin(GL_POINTS);

				if (x>x_min&&x<x_max&&y>y_min&&y<y_max)
				{
					float x_ = x + ICONSIZE / 2;
					float y_ = y + ICONSIZE / 2;
					glVertex2f(x_, y_);
					//count++;
				}
				glEnd();
			}
		}
	}
	
	

}

void BubbleGameEngine::drawBubbleCurve()
{
	//boundingBox();
	int count = 0;
	int x_icon, y_icon;
	int x_min, y_min, x_max, y_max;
	int bubble_size = bubble_clusters_.size();
	for (int i = 0; i < bubble_size;i++)
	{
		bubble_clusters_[i]->getBoundingBox(&x_min, &y_min, &x_max, &y_max);
		for (int x = x_min; x < x_max; x++)
		{
			for (int y = y_min; y < y_max; y++)
			{
				vector<Icon*> cluster_i = bubble_clusters_[i]->return_icon();
				int cluster_i_length = cluster_i.size();
				float* point_to_icon_distance = new float[cluster_i_length];
				float h0_i = hmax;
				for (int k = 0; k < cluster_i_length; k++)
				{
					cluster_i[k]->getIconPosition(&x_icon, &y_icon);
					point_to_icon_distance[k] = radialQuadraticFunction(x, x_icon + CANVASW,y,y_icon + CANVASH);
					h0_i = h0_i + point_to_icon_distance[k];
				}
				float* point_to_icon_weight = new float[cluster_i_length];
				float fx = 0;
				for (int l = 0; l < cluster_i_length; l++)
				{
					point_to_icon_weight[l] = hmax / h0_i;
					fx = fx + (point_to_icon_weight[l]) * (point_to_icon_distance[l]);
				}
				delete[] point_to_icon_distance;
				delete[] point_to_icon_weight;
				glPointSize(5.0f);
				glBegin(GL_POINTS);
				if (0.97 < fx && fx < 1.03)
				{
					glVertex2f(x, y);
					count++;
				}			
				glEnd();
			}
		}
	}
}



/*void BubbleGameEngine::boundingBox()
{	
	int *boungind_box = new int[4];
	//float bounding_x1, float bounding_x2, float bounding_y1, float bounding_y2;
	int bubble_size = bubble_clusters_.size();
	for (int i = 0; i < bubble_size; i++)
	{
		vector<Icon*> cluster_i = (bubble_clusters_[i])->return_icon();
		int cluster_length = cluster_i.size();
		int x_min = cluster_i[0]->returnPosition()[0]; 
		int y_min = cluster_i[0]->returnPosition()[1];
		int x_max = cluster_i[0]->returnPosition()[0];
		int y_max = cluster_i[0]->returnPosition()[1];
		for (int j = 0; j < cluster_length; j++) // for each icon j in cluster i
		{
			int x = cluster_i[j]->returnPosition()[0];
			int y = cluster_i[j]->returnPosition()[1];
			if (x >= x_max)
			{
				x_max = x;
			}
			if (y >= y_max)
			{
				y_max = y;
			}
			if (x <= x_min)
			{
				x_min = x;
			}
			if (y <= y_min)
			{
				y_min = y;
			}
		}
		if (x_min - CLUSTERBOUND < 0) x_min = 0;
		if (y_min - CLUSTERBOUND < 0) y_min = 0;
		if (x_max + CLUSTERBOUND > CANVASW) x_max = CANVASW;
		if (y_max + CLUSTERBOUND > CANVASH) y_max = CANVASH;
		bubble_clusters_[i]->setBoundingBox(x_min - CLUSTERBOUND, y_min - CLUSTERBOUND, x_max + CLUSTERBOUND, y_max + CLUSTERBOUND);

	}
}*/

void BubbleGameEngine::movedDrawCurve(int x, int y)
{
	int count = 0;
	int x_icon, y_icon;
	int xmin, ymin, xmax, ymax;
	//int x, y;  // mouse's coordinates
	vector<BubbleCluster*>::iterator cluster_iter;
	for (cluster_iter = bubble_clusters_.begin(); cluster_iter != bubble_clusters_.end(); ++ cluster_iter)
	{
		if ((*cluster_iter)->clusterSlected(x, y))
		{
			(*cluster_iter)->getBoundingBox(&xmin, &ymin, &xmax, &ymax);
			for (int x_box = xmin; x < xmax; x++)
			{
				for (int y_box = ymin; y < ymax; y++)
				{
					if ((*cluster_iter)->clusterSlected(x_box, y_box)) continue;
					vector<Icon*> cluster_i = (*cluster_iter)->return_icon();
					int cluster_i_length = cluster_i.size();
					float* point_to_icon_distance = new float[cluster_i_length];
					float h0_i = hmax;
					for (int k = 0; k < cluster_i_length; k++)
					{
						cluster_i[k]->getIconPosition(&x_icon, &y_icon);
						point_to_icon_distance[k] = radialQuadraticFunction(x, x_icon + CANVASW, y, y_icon + CANVASH);
						h0_i = h0_i + point_to_icon_distance[k];
					}
					float* point_to_icon_weight = new float[cluster_i_length];
					float fx = 0;
					for (int l = 0; l < cluster_i_length; l++)
					{
						point_to_icon_weight[l] = hmax / h0_i;
						fx = fx + (point_to_icon_weight[l]) * (point_to_icon_distance[l]);
					}
					glPointSize(2);
					glBegin(GL_POINTS);
					if (0.97 < fx && fx < 1.03)
					{
						glVertex2f(x, y);
						count++;
					}
					glEnd();
				}
			}

		}
	}
}



void BubbleGameEngine::moveBubble(int x, int y)
{
	float clusterFx;
	bool clusterSelected = false;
	bool iconSelected =  false;
	int x_icon, y_icon;
	vector<BubbleCluster*>::iterator cluster_iter;
	vector<Icon*>::iterator icon_iter;
	for (cluster_iter = bubble_clusters_.begin(); cluster_iter != bubble_clusters_.end(); ++cluster_iter)
	{
		(*cluster_iter)->clusterWeight(x, y);
		(*cluster_iter)->getClusterFx(&clusterFx);
		vector<Icon*> icons_ = (*cluster_iter)->return_icon();
		for (icon_iter = icons_.begin(); icon_iter != icons_.end(); ++icon_iter)
		{
			(*icon_iter)->getIconPosition(&x_icon, &y_icon);
			if ((x > x_icon && x < x_icon + ICONSIZE && y > y_icon && y < y_icon + ICONSIZE) && (iconSelected = false))
			{
				iconSelected = true;
			}
		}
		iconSelected = false;
		if (clusterFx > 1 && iconSelected == false)
		{
			for (icon_iter = icons_.begin(); icon_iter != icons_.end(); ++icon_iter)
			{
				(*icon_iter)->setIconSlected(true);
			}
		}
		/*	clusterSelected = true;
		}
		if (clusterSelected)
		{
			for (icon_iter = ((*cluster_iter)->return_icon()).begin(); icon_iter != ((*cluster_iter)->return_icon()).end(); ++icon_iter)
			{
				(*icon_iter)->getIconPosition(&x_icon, &y_icon);
				x_icon = x_icon + x_move;
				y_icon = y_icon + y_move;
			}
			(*cluster_iter)->findIcon(x, y, x_move, y_move);
		}
		(*cluster_iter)->clearClusterSelected();*/
	}
}

void BubbleGameEngine::findBubble(int x, int y)
{
	for (int i = 0; i < bubble_clusters_.size(); i++)
	{
		bubble_clusters_[i]->flag_ = 0;
		vector<Icon*> icons = bubble_clusters_[i]->return_icon();
		for (int j = 0; j < icons.size(); j++)
		{
			if (icons[j]->iconSel(x, y))
			{
				bubble_clusters_[i]->flag_ = 1;
				//cout << "cluster " << i << " flag is " << bubble_clusters_[i]->flag_ << endl;
			}
		}
		if (bubble_clusters_[i]->flag_ == 1 && bubble_clusters_[i]->return_icon().size() >=3)
		{
			bubble_clusters_[i]->triangulateT();
			//cout << "hehhh" << endl;
			bubble_clusters_[i]->spreading();
			//cout << "hhhhhhhh" << endl;
		}
		if (bubble_clusters_[i]->flag_ == 1 && bubble_clusters_[i]->return_icon().size() == 2)
		{
			bubble_clusters_[i]->spreadingTwo();
			//cout << "hhehwoho" << endl;
		}
	}
}




float radialQuadraticFunction(int x1, int x2, int y1, int y2 )
{
	int r0 = RADIUS;
	// int r1 = r0 * 2;
	float a, b, c; // f(x) = ax^2 + bx + c
	a = 1.0 / (r0*r0);
	b = -4.0 / r0;
	c = 4.0;
	float x = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
	float f_x;
	if (x > 2 * r0) f_x = 0;
	else f_x = a*x*x + b*x + c;
	return f_x;
}
