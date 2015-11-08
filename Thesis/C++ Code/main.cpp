// C/C++ Libraries
#include <iostream>
#include<typeinfo>

#include <Eigen/Eigenvalues>

// Qt Libraries
//#include <QCoreApplication>

// PCL Libraries
#include <pcl/visualization/cloud_viewer.h>
#include <pcl/io/io.h>
#include <pcl/io/pcd_io.h>

#include <pcl/search/search.h>
#include <pcl/search/kdtree.h>
#include <pcl/features/normal_3d.h>

// Local Libraries
#include "regiongrowingsegmentation.h"

//C++ file writing
#include <fstream>

// OpenMP - MultiThreading
#include <omp.h>

#include <math.h>

#include <time.h>

time_t start, end;

int user_data;

void
viewerOneOff(pcl::visualization::PCLVisualizer& viewer)
{
	viewer.setBackgroundColor(0.0, 0.0, 0.0);
	pcl::PointXYZ o;
	o.x = 1.0;
	o.y = 0;
	o.z = 0;
	viewer.addSphere(o, 0.25, "sphere", 0);
	std::cout << "i only run once" << std::endl;
}

void
viewerPsycho(pcl::visualization::PCLVisualizer& viewer)
{
	static unsigned count = 0;
	std::stringstream ss;
	ss << "Once per viewer loop: " << count++;
	viewer.removeShape("text", 0);
	viewer.addText(ss.str(), 200, 300, "text", 0);

	//FIXME: possible race condition here:
	user_data++;
}


int main(int argc, char *argv[])
{

	time(&start);
	//    QCoreApplication a(argc, argv);

	// ------------------------------------------------------------
	// load the point cloud
	// ------------------------------------------------------------
	//    pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZRGBA>);
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
	if (cloud == 0)
	{
		cout << "Error" << endl;
	}
	else
	{
		//Importing Point Cloud
		pcl::io::loadPCDFile("C:/Users/jayre/Documents/Visual Studio 2013/Projects/ConsoleApplication1/Clouds/Snape_0_01.pcd", *cloud);

		//Getting centroid of cloud
		Eigen::Vector4f centroid;

		pcl::compute3DCentroid(*cloud, centroid);

		pcl::PointXYZRGB bottomZ, topZ;

		// Get min and max 3D point for cloud
		pcl::getMinMax3D(*cloud, bottomZ, topZ);

		cout << "File loaded" << endl;
		cout << "Point count " << cloud->points.size() << endl;
		std::cout << "Centroid: " << std::endl << centroid << std::endl;
		std::cout << "Minimum Z height is: " << bottomZ.z << std::endl;
		std::cout << "Maximum Z height is: " << topZ.z << std::endl;
	}


	// ------------------------------------------------------------
	// segment the point cloud
	// ------------------------------------------------------------
	P::RegionGrowingSegmentation region_segmentor;

	int knn = 10;
	float smoothness_threshold = 5.0;
	float curvature_threshold = 5.0;
	int min_cluster_size = 500;
	int max_cluster_size = 1e8;

	region_segmentor.setParameters(knn, smoothness_threshold, curvature_threshold, min_cluster_size, max_cluster_size);
	region_segmentor.segment(cloud);

	//Timing
	time(&end);//Stops timing
	double dif = difftime(end, start);
	double min = (dif / 60);
	int rMin = min;
	double sec = (min - rMin) * 60;
	cout << "Elapsed time is " << rMin << " minutes " << sec << " seconds" << endl;

	// ------------------------------------------------------------
	// view the point cloud
	// ------------------------------------------------------------
	pcl::visualization::CloudViewer viewer("Cloud Viewer");

	//blocks until the cloud is actually rendered
	viewer.showCloud(cloud);

	//use the following functions to get access to the underlying more advanced/powerful
	//PCLVisualizer

	// This will only get called once
	//viewer.runOnVisualizationThreadOnce (viewerOneOff);

	// This will get called once per visualization iteration
	//viewer.runOnVisualizationThread (viewerPsycho);


	while (!viewer.wasStopped())
	{
		//you can also do cool processing here
		//FIXME: Note that this is running in a separate thread from viewerPsycho
		//and you should guard against race conditions yourself...
		user_data++;
	}

	//    return a.exec();
}
