#pragma once

#include <math.h>
#include <iostream>
#include <vector>
//#include "triangle.h"
#include "common.h"
#include "Vector.h"
#include <Eigen/dense>

using namespace std;
using namespace Eigen;

#define CANVASW 800		// drawing canvas bound width
#define CANVASH 600	    // height
#define ICONSIZE 25    	// icon size
#define RADIUS  30		// based bubble radius
#define SMALL_THREASHOLD 55 // small theashold for clustering
#define LARGE_THREASHOLD 70 // large theashold for spliting
#define hmax 4 
#define CLUSTERBOUND 50 // cluster bounding box
#define K 10 // icon number
#define NEGATIVEAFFECT 0.3 // negative affect from other clusters
#define SPREAD_DISTANCE 15
#define SPREAD_DISTANCE_TWO 10

class Point2d
{
public:
	Point2d() {}
	Point2d(int x, int y) { 
		this->x_ = x; 
		this->y_ = y;
	}
	~Point2d() {}

public:
	int x_;
	int y_;
public:
	//int x_;
	//int y_;
	int index;
	int flag_;
};


class Icon
{
public:
	Icon(float x, float y, float w, float h) : x_(x), y_(y), w_(w), h_(h) {}	// x, y::position, w, h::size
	~Icon();

public:
	void draw();
	void moveIcon(int x, int y, int x_move, int y_move);
	void getIconPosition(int *x, int *y)
	{
		*x = x_;
		*y = y_;
	}

	void setGroupLable(int groupLable)
	{
		group_label_ = groupLable;
	}
	void getGroupLable(int *groupLable)
	{
		*groupLable = group_label_;
	}

	void clearSlected()
	{
		iconSlected = false;
	}
	void setIconSlected(bool icon_selected)
	{
		iconSlected = icon_selected;
	}
	void setIconWeight(float iconWeight)
	{
		icon_weight_ = iconWeight;
	}
	void getIconWeight(float *iconWeight)
	{
		*iconWeight = icon_weight_;
	}
	bool iconSel(int x, int y)
	{
		if ((x_ < x) && x < (x_ + ICONSIZE) && (y_ < y) && (y < y_ + ICONSIZE))
		{
			return true;
		}
		else
			return false;
	}
	void setIconPosition(int x, int y)
	{
		this->x_ = x;
		this->y_ = y;
	}

private:
	bool iconSlected;
	float x_, y_;		// icon position
	float w_, h_;		// icon size
	int group_label_;	// which group
	float icon_weight_;
};




class BubbleCluster
{
public:
	BubbleCluster(vector<Icon*> icons); // vetors are sequence containers arrays that
	BubbleCluster(Icon* icon);          // can change in size
	BubbleCluster();
	~BubbleCluster();

public:
	void draw();
	void addIcon(Icon* icon)
	{
		this->icons_.push_back(icon);
	}
	void addIcons(vector<Icon*> newicons)
	{
		vector<Icon*>::iterator icon_iter;
		for (icon_iter = newicons.begin(); icon_iter != newicons.end(); ++icon_iter)
		{
			this->icons_.push_back(*icon_iter);
		}  // add element at the end 
	}
	void removeIcon(Icon* icon)
	{
		vector<Icon*>::iterator icon_iter;
		// begin(): return iterator to beginning, similar to end()
		for (icon_iter = icons_.begin(); icon_iter != icons_.end(); ++icon_iter)
		{
			if (*icon_iter == icon)
			{
				this->icons_.erase(icon_iter); // erase elements
				return;
			}
		}
	}
	void findIcon(int x, int y, int x_move, int y_move);
	vector<Icon*> return_icon()
	{
		return icons_;
	}
	void clearSlected();
	void boundingBox();
	bool clusterSlected(int x, int y);
	void clusterWeight(int x, int y);
	void BubbleCluster::pointToClusterValue(int x, int y);


	void setGroupLable(int groupLable)
	{
		cluster_group_ = groupLable;
	}
	void getGroupLable(int *groupLable)
	{
		*groupLable = cluster_group_;
	}

	void setClusterFx(float clusterFx)
	{
		cluster_fx_ = clusterFx;
	}
	void getClusterFx(float *clusterFx)
	{
		*clusterFx = cluster_fx_;
	}
	void getBoundingBox(int *x_min, int *y_min, int *x_max, int *y_max)
	{
		*x_min = x_min_;
		*y_min = y_min_;
		*x_max = x_max_;
		*y_max = y_max_;
	}

	void setBoundingBox(int x_min,int y_min, int x_max, int y_max)
	{
		x_min_ = x_min;
		y_min_ = y_min;
		x_max_ = x_max;
		y_max_ = y_max;
	}
	void clearClusterSelected()
	{
		clear_cluster_selected = false;
	}



	void triangulateInit(triangulateio * io);
	void triangulateT();
	void trianguleReport(struct triangulateio *io, int markers, int reporttriangles, int reportneighbors, int reportsegments,
		int reportedges, int reportnorms);

	
	void spreadingTwo();
	void spreading();


public:
	int flag_;

private:
	vector<Icon*> icons_;
	int cluster_group_;
	float cluster_fx_;
	int x_min_, y_min_, x_max_, y_max_;
	bool clear_cluster_selected;
	triangulateio out_tri_;
};


class BubbleGameEngine
{
public:
	BubbleGameEngine();
	~BubbleGameEngine();
	

public:
	void init();	// initialize, randomly generate a set of bubble clusters
	void draw();
	void findCluster(int x, int y, int x_move, int y_move);
	//int * returnIconPosition();
	//int calculateDistance(int x, int y, int x_, int y_);
	void removeCluster(vector<BubbleCluster*>* bclusters, BubbleCluster* cluster)
	{
		vector<BubbleCluster*>::iterator cluster_iter;
		// begin(): return iterator to beginning, similar to end()
		for (cluster_iter = bclusters->begin(); cluster_iter != bclusters->end(); ++cluster_iter)
		{
			if (*cluster_iter == cluster)
			{
				bclusters->erase(cluster_iter); // erase elements
				return;
			}
		}
	}
	void removeCluster(BubbleCluster* cluster)
	{
		vector<BubbleCluster*>::iterator cluster_iter;
		// begin(): return iterator to beginning, similar to end()
		for (cluster_iter = bubble_clusters_.begin(); cluster_iter != bubble_clusters_.end(); ++cluster_iter)
		{
			if (*cluster_iter == cluster)
			{
				bubble_clusters_.erase(cluster_iter); // erase elements
				return;
			}
		}
	}
	
	bool ifClusterMerge(vector<Icon*> icons_i, vector<Icon*> icons_j);

	vector<BubbleCluster*> BubbleGameEngine::clusteringInitial(vector<BubbleCluster*> clusters);
	void clustering_spliting();
	void clearSlected();
	
	int calculateDistance(Icon* icon_i, Icon* icon_j);
	void drawBubbleCurve();
	//void boundingBox();
	//float radialQuadraticFunction(int x1, int x2, int y1, int y2);
	void drawCurve();
	void drawCurve2();
	void drawBoundingBox();
	void movedDrawCurve(int x, int y);
	void moveBubble(int x, int y);
	
	void findBubble(int x, int y);

private:
	vector<BubbleCluster*> bubble_clusters_;
};

float radialQuadraticFunction(int x1, int x2, int y1, int y2);