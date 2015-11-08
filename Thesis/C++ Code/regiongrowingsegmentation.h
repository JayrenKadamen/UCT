#ifndef REGIONGROWINGSEGMENTATION_H
#define REGIONGROWINGSEGMENTATION_H

// C/C++ Libraries
#include <iostream>
//For deleting files
#include <stdio.h>

// PCL Libraries
#include <pcl/visualization/cloud_viewer.h>
//PCL Segmentation
#include <pcl/features/normal_3d.h>
#include <pcl/io/pcd_io.h>
#include <pcl/segmentation/region_growing.h>
//PCL Min/Max X,Y,Z
#include <pcl/point_types.h>
#include <pcl/common/common.h>

//Math library
#include <math.h>
#define PI 3.14159265

//EigenValue Stuff
#include <pcl/common/pca.h>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

//C++ file writing
#include <fstream>

// OpenMP - MultiThreading
#include <omp.h>

//To delete files in a folder
#include <boost\filesystem.hpp>

namespace P
{

	struct SegmentationMetrics
	{
		double sum; // see getCurvatureHomogeneityMeasure
		double seg_stddev; // see getCurvatureHomogeneityMeasure
		double NPSS;
		double IntraRegionUniformity;
		double AverageAbsoluteDeviation;
		double CurvatureHeterogeneity;
		double P;
		double RGBHeterogeneity;
		double RGBIntraRegionUniformity;
		double RGBAverageAbsoluteDeviation;
	};


	struct MetricsVariables
	{
		int segcount;
		int point_num;

		double feature_average;
		double Max_f;
		double Min_f;
		double df;
		double N;
		double N_abs;

		double Max_frgb;
		double Min_frgb;
		double dfrgb;
		double Nrgb;
		double N_devrgb;
	};

	/** This class performs a region growing segmentation
	*/
	class RegionGrowingSegmentation
	{
	public:
		RegionGrowingSegmentation();
		~RegionGrowingSegmentation();


		// type definitions
		typedef pcl::PointCloud<pcl::PointXYZRGBA> ColorACloud;
		typedef pcl::PointCloud<pcl::PointXYZRGB> ColorCloud;
		typedef pcl::PointCloud<pcl::PointXYZI> GrayCloud;
		typedef pcl::PointCloud<pcl::PointXYZ> MonochromeCloud;

		typedef pcl::PointCloud <pcl::Normal> NormalCloud;


		/** \brief Segment a cloud.
		* \param[in] cloud RGB point cloud
		*/
		void segment(const ColorCloud::ConstPtr &cloud);

		/** \brief Show a cloud.
		* \param[in] cloud RGB point cloud
		*/
		void segment(const ColorACloud::ConstPtr &cloud);

		/** \brief Show a cloud.
		* \param[in] cloud XYZI point cloud
		*/
		void segment(const GrayCloud::ConstPtr &cloud);


		/** \brief Show a cloud.
		* \param[in] cloud XYZ point cloud
		*/
		void segment(const MonochromeCloud::ConstPtr &cloud);

		/** \brief Set parameters
		* \param[in] to be explained
		*/
		void setParameters(const int& knn,
			const float& smoothness_threshold, const float& curvature_threshold,
			const int& min_cluster_size, const int& max_cluster_size)
		{
			m_knn = knn;
			m_smoothness_threshold = smoothness_threshold;
			m_curvature_threshold = curvature_threshold;

			m_min_cluster = min_cluster_size;
			m_max_cluster = max_cluster_size;
		}

	public: Eigen::Vector4f centroid;

	protected:
		double background_average(const pcl::PointCloud<pcl::Normal>::ConstPtr &cloud, int num_p);

		double region_average(pcl::PointCloud<pcl::Normal>::Ptr &cloudn, int num_p);

		double compute_heterogeneity(double b_average, double r_average);

		double min(pcl::PointCloud<pcl::Normal>::Ptr &merged, int num_p);

		double max(pcl::PointCloud<pcl::Normal>::Ptr &merged, int num_p);

		double compute_average(int num_p, pcl::PointCloud<pcl::Normal>::Ptr &cloud);

		double compute_v(double mean, pcl::PointCloud<pcl::Normal>::Ptr &cloud);

		double min_rgb(const pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr &merged, int num_p);

		double max_rgb(const pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr &merged, int num_p);

		double compute_average_rgb(int num_p, pcl::PointCloud<pcl::PointXYZRGB>::Ptr &cloud);

		double compute_v_rgb(double mean, int num_p, pcl::PointCloud<pcl::PointXYZRGB>::Ptr &cloud);

		double compute_abs(double mean, int num_p, pcl::PointCloud<pcl::Normal>::Ptr &cloud);

		double compute_absrgb(double mean, int num_p, pcl::PointCloud<pcl::PointXYZRGB>::Ptr &cloud);

		double background_averagergb(const pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr &cloudrgb, int num_p);

		double region_averagergb(const pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr &cloudrgb, int num_p);

		double compute_heterogeneityrgb(double b_averagergb, double r_averagergb);


		// --------------------------------------------
		// calculate the normals for the point cloud and store the normals in normalcloud
		template <class Cloud, class PointType, class Tree>
		void getNormals(Cloud &cloud, Tree &tree, NormalCloud::Ptr &normalcloud)
		{
			std::cout << "--------------------------------------------" << std::endl;
			std::cout << "Computing Normals ..." << std::endl;

			pcl::NormalEstimation<PointType, pcl::Normal> normal_estimator;
			normal_estimator.setSearchMethod(tree);
			normal_estimator.setInputCloud(cloud);
			normal_estimator.setKSearch(m_knn);
			normal_estimator.compute(*normalcloud);


			//pcl::visualization::PCLVisualizer viewer("Cloud and Normals");

			//pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(cloud);
			//viewer.addPointCloud<pcl::PointXYZRGB>(cloud, rgb, "sample cloud");

			//viewer.setBackgroundColor(0.0, 0.0, 0.0);
			//viewer.addPointCloudNormals<pcl::PointXYZRGB, pcl::Normal>(cloud, normalcloud, 400, 0.5, "normals");

			//while (!viewer.wasStopped())
			//{
			//	viewer.spinOnce();
			//}


			std::cout << "Normals Computed ... " << normalcloud->points.size() << std::endl << std::endl;
		}

		

		// --------------------------------------------
		// perform the region growing segmentation
		template <class Cloud, class PointType, class Tree, class Indeces>
		void getSegments(Cloud &cloud, Tree &tree, NormalCloud::Ptr &normalcloud, int knn, Indeces &clusters)
		{
			std::cout << "--------------------------------------------" << std::endl;
			std::cout << "Extracting Segments ..." << std::endl;

			float smoothness_threshold = float(m_smoothness_threshold / 180.0 * M_PI);
			float curvature_threshold = float(m_curvature_threshold / 180.0 * M_PI);

			pcl::RegionGrowing<pcl::PointXYZRGB, pcl::Normal> reg;
			reg.setMinClusterSize(m_min_cluster);
			reg.setMaxClusterSize(m_max_cluster);
			reg.setSearchMethod(tree);
			reg.setNumberOfNeighbours(knn);
			reg.setInputCloud(cloud);
			//reg.setIndices (indices);
			reg.setInputNormals(normalcloud);
			reg.setSmoothnessThreshold(smoothness_threshold);
			reg.setCurvatureThreshold(curvature_threshold);

			//        Indeces clusters;
			reg.extract(clusters);

			pcl::PointCloud <pcl::PointXYZRGB>::Ptr colored_cloud = reg.getColoredCloud();
			pcl::PCDWriter writer;
			writer.write("colouredcloud.pcd", *colored_cloud, false);

			std::cout << "Segments Extracted ... " << clusters.size() << std::endl << std::endl;
		}

		// --------------------------------------------
		// write Segments to separate files
		template <class Cloud, class PointType, class Indeces>
		void writeSegmentsToSeparateFiles(Cloud &cloud, NormalCloud::Ptr &normalcloud, Indeces &clusters)
		{

			// declare file writers
			pcl::PCDWriter writer;

			std::cout << "Writing " << clusters.size() << " segments to separate files ..." << clusters.size() << std::endl;
			int segindex = 0;

			// Delete old files
			//--------------------------------------------------------------------------
			// CSV Files
			if (remove("minmax/MinMaxZ.csv") != 0)
				std::cout << "No file to delete - Continuing" << std::endl;
			else
				std::cout << "Old CSV files have been successfully cleared" << std::endl;

			//if (remove("eValues/minEigenvalue.csv") != 0)
			//	std::cout << "No file to delete - Continuing" << std::endl;
			//else
			//	std::cout << "Old CSV files have been successfully cleared" << std::endl;

			if (remove("HorizontalPlanes/H_Plane_Data.csv") != 0)
				std::cout << "No file to delete - Continuing" << std::endl;
			else
				std::cout << "Old CSV files have been successfully cleared" << std::endl;

			if (remove("VerticalPlanes/V_Plane_Data.csv") != 0)
				std::cout << "No file to delete - Continuing" << std::endl;
			else
				std::cout << "Old CSV files have been successfully cleared" << std::endl;

			if (remove("OtherSegments/Other_Segment_Data.csv") != 0)
				std::cout << "No file to delete - Continuing" << std::endl;
			else
				std::cout << "Old CSV files have been successfully cleared" << std::endl;

			// PCD Files
			// Remove all files in folder, sadly also the folder :(
			if (boost::filesystem::remove_all("HorizontalPlanes/") != 0)
				std::cout << "No file to delete - Continuing" << std::endl;
			else
				std::cout << "HorizontalPlanes directory has been successfully deleted" << std::endl;
			// Need to recreate the folder that was just deleted
			boost::filesystem::create_directory("HorizontalPlanes");

			if (boost::filesystem::remove_all("VerticalPlanes/") != 0)
				std::cout << "No file to delete - Continuing" << std::endl;
			else
				std::cout << "VerticalPlanes directory has been successfully deleted" << std::endl;
			boost::filesystem::create_directory("VerticalPlanes");

			if (boost::filesystem::remove_all("OtherSegments/") != 0)
				std::cout << "No file to delete - Continuing" << std::endl;
			else
				std::cout << "OtherSegments directory has been successfully deleted" << std::endl;
			boost::filesystem::create_directory("OtherSegments");

			if (boost::filesystem::remove_all("SegmentData/") != 0)
				std::cout << "No file to delete - Continuing" << std::endl;
			else
				std::cout << "SegmentData directory has been successfully deleted" << std::endl;
			boost::filesystem::create_directory("SegmentData");

			if (boost::filesystem::remove_all("normals/") != 0)
				std::cout << "No file to delete - Continuing" << std::endl;
			else
				std::cout << "normals directory has been successfully deleted" << std::endl;
			boost::filesystem::create_directory("normals");

			if (boost::filesystem::remove_all("segments/") != 0)
				std::cout << "No file to delete - Continuing" << std::endl;
			else
				std::cout << "segments directory has been successfully deleted" << std::endl;
			boost::filesystem::create_directory("segments");

			// Creates 3D points min_pt and max_pt
			pcl::PointXYZRGB min_pt1, max_pt1;
			pcl::PointXYZRGB bottomZ, topZ;


			//Loops through all segments to find range of z values
			for (std::vector<pcl::PointIndices>::const_iterator it = clusters.begin(); it != clusters.end(); ++it)
			{
				pcl::PointCloud<PointType>::Ptr point_cloud_seg1(new pcl::PointCloud<PointType>);

				for (std::vector<int>::const_iterator pit = it->indices.begin(); pit != it->indices.end(); pit++)
				{	//Pushing each segment to the allocated variable in memory
					point_cloud_seg1->points.push_back(cloud->points[*pit]);
				}

				// setup the point cloud
				point_cloud_seg1->width = point_cloud_seg1->points.size();
				point_cloud_seg1->height = 1;
				point_cloud_seg1->is_dense = true;

				// Get min and max 3D point for cloud
				pcl::getMinMax3D(*point_cloud_seg1, min_pt1, max_pt1);

				if (min_pt1.z<bottomZ.z)
				{
					bottomZ.z = min_pt1.z;
				}

				if (max_pt1.z > topZ.z)
				{
					topZ.z = max_pt1.z;
				}
			}
			std::cout << "Height ranges from: " << bottomZ.z << " to " << topZ.z << std::endl;

			for (std::vector<pcl::PointIndices>::const_iterator it = clusters.begin(); it != clusters.end(); ++it)
			{
				pcl::PointCloud<pcl::Normal>::Ptr normal_cloud_seg(new pcl::PointCloud<pcl::Normal>);
				pcl::PointCloud<PointType>::Ptr point_cloud_seg(new pcl::PointCloud<PointType>);

				//        pcl::PointCloud<MyPointType>::Ptr first_seg (new pcl::PointCloud<MyPointType>);

				for (std::vector<int>::const_iterator pit = it->indices.begin(); pit != it->indices.end(); pit++)
				{
					normal_cloud_seg->points.push_back(normalcloud->points[*pit]);
					point_cloud_seg->points.push_back(cloud->points[*pit]);
				}

				// setup the normal cloud
				normal_cloud_seg->width = normal_cloud_seg->points.size();
				normal_cloud_seg->height = 1;
				normal_cloud_seg->is_dense = true;

				// setup the point cloud
				point_cloud_seg->width = point_cloud_seg->points.size();
				point_cloud_seg->height = 1;
				point_cloud_seg->is_dense = true;

				//m_mv.segcount = clusters.size();
				//m_mv.point_num = countOfPointsInAllSegments<pcl::PointXYZRGB>(m_mv.segcount);
				//cout << "Number of points in all segments: " << m_mv.point_num << endl;

				// Creates 3D points min_pt and max_pt
				pcl::PointXYZRGB min_pt, max_pt;

				// ------------------------------------------------------------------------------------------------------
				// My PCA
				// ------------------------------------------------------------------------------------------------------
				pcl::PCA<pcl::PointXYZRGB> pca;
				pca.setInputCloud(point_cloud_seg);

				Eigen::MatrixXf e_vec = pca.getEigenVectors();
				Eigen::MatrixXf e_val = pca.getEigenValues();
				Eigen::MatrixXf eMin = pca.getEigenValues();

				Eigen::MatrixXf::Index minRow, minCol;
				float min = eMin.minCoeff(&minRow, &minCol);

				std::cout << "PCA EigenValues: " << std::endl << e_val << std::endl << std::endl;
				std::cout << "PCA EigenVectors: " << std::endl << e_vec << std::endl << std::endl;
				std::cout << "Minimum Eigenvalue: " << min << " at " << minCol << minRow << std::endl << std::endl;
				std::cout << "Corresponding EigenVector: " << std::endl << e_vec.col(minRow) << std::endl << std::endl;

				Eigen::Vector3f final_e_vec = e_vec.col(minRow);
				std::cout << "FINAL E-VEC is: " << final_e_vec << std::endl;
				//Angle between normal and vertical
				std::cout << "Angle is: " << acos(final_e_vec(2)) * 180.0 / PI << " degrees" << std::endl;
				double angle = acos(final_e_vec(2)) * 180.0 / PI;

				// Get min and max 3D point for cloud
				pcl::getMinMax3D(*point_cloud_seg, min_pt, max_pt);

				//Ensuring the ratio of smallest eigenvalue to remaining values is good - i.e. the normal to the plane is well defined
				for (int index = 0; index < 3; index = index + 1)
				{
					if (index != minRow)
					{
						float ratio = e_val(minRow) / e_val(index);
						std::cout << "Eigenvalue: " << index << ": " << e_val(index) << std::endl;
						std::cout << "Ratio: " << index << ": " << ratio << std::endl;

						if (ratio >= 1)
						{
							std::cout << "WARNING: Poor Ratio - Normal to the plane is ill defined" << std::endl;
						}
						else{}
					}
					else{}
				}


				// Horizontal Segment
				if (angle > 0.0 && angle < 20.0) {
					//Write Horizontal Segment to PCD file
					stringstream filename_horizontal;
					filename_horizontal << "HorizontalPlanes/H_segment_" << segindex << ".pcd";
					writer.write<PointType>(filename_horizontal.str(), *point_cloud_seg, false);

					float average_pt = (((min_pt.z + abs(bottomZ.z)) + (max_pt.z + abs(bottomZ.z))) / 2);
					std::cout << "Average Z for segment: " << average_pt << std::endl;
					float min_point = (min_pt.z + abs(bottomZ.z));
					float max_point = (max_pt.z + abs(bottomZ.z));

					if (average_pt > -1 && average_pt < 0.15)
					{

						//Create new csv file
						std::ofstream outfile;
						outfile.open("SegmentData/H_Plane_Data.csv", std::ios_base::app);
						//Write to file
						outfile << "H_Seg_" << segindex <<
							",PtCount," << (*point_cloud_seg).width <<
							",Angle," << angle <<
							",Min_pt," << min_point <<
							",Max_pt," << max_point <<
							",Average Height," << average_pt <<
							",Floor" <<
							std::endl;
					}

					else if (average_pt > 0.45 && average_pt < 0.53)
					{
						//Create new csv file
						std::ofstream outfile;
						outfile.open("SegmentData/H_Plane_Data.csv", std::ios_base::app);
						//Write to file
						outfile << "H_Seg_" << segindex <<
							",PtCount," << (*point_cloud_seg).width <<
							",Angle," << angle <<
							",Min_pt," << min_point <<
							",Max_pt," << max_point <<
							",Average Height," << average_pt <<
							",Chair" <<
							std::endl;
					}

					else if (average_pt > 0.72 && average_pt < 0.80)
					{
						//Create new csv file
						std::ofstream outfile;
						outfile.open("SegmentData/H_Plane_Data.csv", std::ios_base::app);
						//Write to file
						outfile << "H_Seg_" << segindex <<
							",PtCount," << (*point_cloud_seg).width <<
							",Angle," << angle <<
							",Min_pt," << min_point <<
							",Max_pt," << max_point <<
							",Average Height," << average_pt <<
							",Desk" <<
							std::endl;
					}

					else if (average_pt > 2.5 && point_cloud_seg->size() > 10000)
					{
						//Create new csv file
						std::ofstream outfile;
						outfile.open("SegmentData/H_Plane_Data.csv", std::ios_base::app);
						//Write to file
						outfile << "H_Seg_" << segindex <<
							",PtCount," << (*point_cloud_seg).width <<
							",Angle," << angle <<
							",Min_pt," << min_point <<
							",Max_pt," << max_point <<
							",Average Height," << average_pt <<
							",Ceiling" <<
							std::endl;
					}

					else
					{
						//Create new csv file
						std::ofstream outfile;
						outfile.open("SegmentData/H_Plane_Data.csv", std::ios_base::app);
						//Write to file
						outfile << "H_Seg_" << segindex <<
							",PtCount," << (*point_cloud_seg).width <<
							",Angle," << angle <<
							",Min_pt," << min_point <<
							",Max_pt," << max_point <<
							",Average Height," << average_pt <<
							",Unknown" << 
							std::endl;
					}

				}

				else if (angle > 160.0 && angle < 180.0){
					//Write Horizontal Point Cloud to PCD file
					stringstream filename_horizontal;
					filename_horizontal << "HorizontalPlanes/H_segment_" << segindex << ".pcd";
					writer.write<PointType>(filename_horizontal.str(), *point_cloud_seg, false);

					float average_pt = (((min_pt.z + abs(bottomZ.z)) + (max_pt.z + abs(bottomZ.z))) / 2);
					std::cout << "Average Z for segment: " << average_pt << std::endl;
					float min_point = (min_pt.z + abs(bottomZ.z));
					float max_point = (max_pt.z + abs(bottomZ.z));

					if (average_pt > -1 && average_pt < 0.15)
					{
						//Create new csv file
						std::ofstream outfile;
						outfile.open("SegmentData/H_Plane_Data.csv", std::ios_base::app);
						//Write to file
						outfile << "H_Seg_" << segindex <<
							",PtCount," << (*point_cloud_seg).width <<
							",Angle," << angle <<
							",Min_pt," << min_point <<
							",Max_pt," << max_point <<
							",Average Height," << average_pt <<
							",Floor" <<
							std::endl;
					}

					else if (average_pt > 0.45 && average_pt < 0.53)
					{
						//Create new csv file
						std::ofstream outfile;
						outfile.open("SegmentData/H_Plane_Data.csv", std::ios_base::app);
						//Write to file
						outfile << "H_Seg_" << segindex <<
							",PtCount," << (*point_cloud_seg).width <<
							",Angle," << angle <<
							",Min_pt," << min_point <<
							",Max_pt," << max_point <<
							",Chair" <<
							std::endl;
					}

					else if (average_pt > 0.70 && average_pt < 0.85)
					{
						//Create new csv file
						std::ofstream outfile;
						outfile.open("SegmentData/H_Plane_Data.csv", std::ios_base::app);
						//Write to file
						outfile << "H_Seg_" << segindex <<
							",PtCount," << (*point_cloud_seg).width <<
							",Angle," << angle <<
							",Min_pt," << min_point <<
							",Max_pt," << max_point <<
							",Desk" <<
							std::endl;
					}

					else
					{
						//Create new csv file
						std::ofstream outfile;
						outfile.open("SegmentData/H_Plane_Data.csv", std::ios_base::app);
						//Write to file
						outfile << "H_Seg_" << segindex <<
							",PtCount," << (*point_cloud_seg).width <<
							",Angle," << angle <<
							",Min_pt," << min_point <<
							",Max_pt," << max_point <<
							",Unknown" << "," << average_pt <<
							std::endl;
					}

				}

				// Vertical Segment
				else if (angle > 70.0 && angle < 110.0){
					//Write Vertical Point Cloud to PCD file
					stringstream filename_vertical;
					filename_vertical << "VerticalPlanes/V_segment_" << segindex << ".pcd";
					writer.write<PointType>(filename_vertical.str(), *point_cloud_seg, false);


					float range = ((max_pt.z + abs(bottomZ.z)) - (min_pt.z + abs(bottomZ.z)));
					std::cout << "Range for this segment: " << range << std::endl;
					float min_point = (min_pt.z + abs(bottomZ.z));
					float max_point = (max_pt.z + abs(bottomZ.z));

					if (range > 1.8 && range < 2.5)
					{
						//Create new csv file
						std::ofstream outfile;
						outfile.open("SegmentData/V_Plane_Data.csv", std::ios_base::app);
						//Write to file
						outfile << "V_Seg_" << segindex <<
							",PtCount," << (*point_cloud_seg).width <<
							",Angle," << angle <<
							",Min_pt," << min_point <<
							",Max_pt," << max_point <<
							",Door" <<
							std::endl;
					}

					else
					{
						//Create new csv file
						std::ofstream outfile;
						outfile.open("SegmentData/V_Plane_Data.csv", std::ios_base::app);
						//Write to file
						outfile << "V_Seg_" << segindex <<
							",PtCount," << (*point_cloud_seg).width <<
							",Angle," << angle <<
							",Min_pt," << min_point <<
							",Max_pt," << max_point <<
							",Other" << "," << range <<
							std::endl;
					}

				}

				// Other Segments
				else {
					//Write Vertical Point Cloud to PCD file
					stringstream filename_other;
					filename_other << "OtherSegments/Segment_" << segindex << ".pcd";
					writer.write<PointType>(filename_other.str(), *point_cloud_seg, false);

					//Create new csv file
					std::ofstream outfile;
					outfile.open("SegmentData/Other_Segment_Data.csv", std::ios_base::app);
					//Write to file
					outfile << "Other_Seg_" << segindex <<
						",Point_Count," << (*point_cloud_seg).width <<
						",Angle" << "," << angle <<
						",Min e-val," << e_val(minRow) <<
						",Min e-vec," << final_e_vec(0) << "," << final_e_vec(1) << "," << final_e_vec(2) <<
						",Min_pt," << min_pt.x << "," << min_pt.y << "," << (min_pt.z + abs(bottomZ.z)) <<
						",Max_pt," << max_pt.x << "," << max_pt.y << "," << (max_pt.z + abs(bottomZ.z)) <<
						std::endl;
				}

				std::cout << "-------------------------------------------------------------" << std::endl;
				//---------------------------------------------------------------------------------------------------------
				//---------------------------------------------------------------------------------------------------------
				// Getting the minimum and maximum z values for each segment
				//---------------------------------------------------------------------------------------------------------

				pcl::PointXYZRGB minPt, maxPt;
				pcl::getMinMax3D(*point_cloud_seg, minPt, maxPt);

				std::cout << "Minimum and Maximum Z for Segment " <<
					segindex << " is: " << minPt.z + abs(bottomZ.z) <<
					" & " << maxPt.z + abs(bottomZ.z) << std::endl;


				//Write
				std::ofstream outfile1;
				outfile1.open("minmax/MinMaxZ.csv", std::ios_base::app);
				outfile1 << "Segment_" << segindex << "," << "MinMax, " << minPt.z + abs(bottomZ.z) << "," << maxPt.z + abs(bottomZ.z) << std::endl;

				/*stringstream filename_minmax;
				filename_minmax << "minmax/minMaxSegment_" << segindex << ".csv";

				myfile.open(filename_minmax.str());
				myfile << "Segment_" << segindex << "," << minPt.z << "," << maxPt.z;
				myfile.close();
				*/
				//---------------------------------------------------------------------------------------------------------

				stringstream filename_normal;
				filename_normal << "normals/segment_" << segindex << ".pcd";
				writer.write<pcl::Normal>(filename_normal.str(), *normal_cloud_seg, false);

				stringstream filename_cloud;
				filename_cloud << "segments/segment_" << segindex << ".pcd";
				writer.write<PointType>(filename_cloud.str(), *point_cloud_seg, false);

				++segindex;

			}
		}

		// --------------------------------------------
		// write Segments to separate files
		template <class PointType>
		void writeSegmentsToSeparateFiles(const int& segcount, const std::string& filenameout)
		{
			pcl::PointCloud<MyPointType> cloud;
			//        if (pcl::io::loadPCDFile<MyPointType> ("rgb_segment_0.pcd", segment) == -1)
			//        {
			//            PCL_ERROR ("Couldn't read file test_pcd.pcd \n");
			//        }

			// FIXME: num_seg should start at 1
			for (int segindex = 0; segindex < segcount; ++segindex)
			{
				stringstream filename_cloudn;

				filename_cloudn << "segments/segment_" << segindex << ".pcd";
				pcl::PointCloud<PointType> cloudn;

				if (pcl::io::loadPCDFile<PointType>(filename_cloudn.str(), cloudn) == -1)
				{
					PCL_ERROR("Couldn't read file test_pcd.pcd \n");
				}
				else
				{
					// create a cloud of MyPointType
					pcl::PointCloud<MyPointType> cloud_seg;
					cloud_seg.width = cloudn.width;
					cloud_seg.height = 1;
					cloud_seg.is_dense = true;
					cloud_seg.points.resize(cloud_seg.width * cloud_seg.height);

					// populate the cloud
					for (size_t i = 0; i < cloudn.points.size(); ++i)
					{
						cloud_seg.points[i].x = cloudn.points[i].x;
						cloud_seg.points[i].y = cloudn.points[i].y;
						cloud_seg.points[i].z = cloudn.points[i].z;
						cloud_seg.points[i].rgb = cloudn.points[i].rgb;
						cloud_seg.points[i].segment_id = segindex;
					}
					cloud += cloud_seg;
				}
			}
			pcl::PCDWriter writer;
			writer.write(filenameout, cloud, false);
		}


		// --------------------------------------------
		// sums the count of points in all the segments
		template <class PointType>
		int countOfPointsInAllSegments(const int& segcount)
		{
			int point_num = 0;
			for (int segindex = 0; segindex < segcount; ++segindex)
			{
				stringstream filename_cloud;
				filename_cloud << "segments/segment_" << segindex << ".pcd";
				pcl::PointCloud<PointType>::Ptr cloud(new pcl::PointCloud<PointType>);
				if (pcl::io::loadPCDFile<PointType>(filename_cloud.str(), *cloud) == -1)
				{
					PCL_ERROR("Couldn't read file test_pcd.pcd \n");
				}
				else
				{
					point_num += cloud->points.size();
				}
			}
			return point_num;
		}

		// preprocessing for metrics
		template <class Cloud>
		void metricsPreProcess(Cloud& cloud, NormalCloud::Ptr& normalcloud)
		{
			// normalised segmentation standard deviation
			// FIXME: check this
			// feature value average for the whole point cloud
			m_mv.feature_average = background_average(normalcloud, normalcloud->points.size());
			m_mv.Max_f = max(normalcloud, normalcloud->points.size());
			m_mv.Min_f = min(normalcloud, normalcloud->points.size());
			m_mv.df = m_mv.Max_f - m_mv.Min_f;
			m_mv.N = (m_mv.point_num * m_mv.df * m_mv.df) / 2.0;
			m_mv.N_abs = (m_mv.point_num* fabs(m_mv.df));

			m_mv.Max_frgb = max_rgb(cloud, cloud->points.size());
			m_mv.Min_frgb = min_rgb(cloud, cloud->points.size());
			m_mv.dfrgb = m_mv.Max_frgb - m_mv.Min_frgb;
			m_mv.Nrgb = (m_mv.point_num * (m_mv.dfrgb * m_mv.dfrgb)) / 2.0;
			m_mv.N_devrgb = (m_mv.point_num * (fabs(m_mv.dfrgb)));
		}

		// get curvature metrics
		void getCurvatureMetrics(NormalCloud::Ptr &normalcloud);

		// --------------------------------------------
		// get colour metrics
		template <class Cloud, class PointType>
		void getColourMetrics(Cloud cloud, NormalCloud::Ptr &normalcloud)
		{
			double sum = 0.0;
			double dev_sumrgb = 0.0;
			for (int segindex = 0; segindex < m_mv.segcount; ++segindex)
			{
				stringstream filename_cloud;

				filename_cloud << "segments/segment_" << segindex << ".pcd";
				//std::cout<<"points on segment_"<<num_seg<<std::endl;

				pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloudrgb(new pcl::PointCloud<pcl::PointXYZRGB>);
				if (pcl::io::loadPCDFile<PointType>(filename_cloud.str(), *cloudrgb) == -1)
				{
					PCL_ERROR("Couldn't read PCD file \n");
				}
				else
				{
					double mean_frgb = compute_average_rgb(cloudrgb->points.size(), cloudrgb);
					double variance = compute_v_rgb(mean_frgb, cloudrgb->points.size(), cloudrgb);
					double abs_devrgb = compute_absrgb(mean_frgb, cloudrgb->points.size(), cloudrgb);

					double max_vrgb = (m_mv.dfrgb * m_mv.dfrgb) / 2.0;
					sum += (cloudrgb->points.size() * (variance) / m_mv.Nrgb);
					dev_sumrgb += (cloudrgb->points.size() * abs_devrgb / m_mv.N_devrgb);
				}
			}

			m_metrics.RGBIntraRegionUniformity = 1 - sum;
			m_metrics.RGBAverageAbsoluteDeviation = 1 - dev_sumrgb;

			cout << "--------------------------------------------" << endl;
			cout << "Colour Metrics" << endl;
			cout << "rgb: intra-Region Uniformity = " << m_metrics.RGBIntraRegionUniformity << endl;
			cout << "rgb: average absolute deviation = " << m_metrics.RGBAverageAbsoluteDeviation << endl;


			double b_average = background_average(normalcloud, normalcloud->points.size());
			double b_averagergb = background_averagergb(cloud, cloud->points.size());

			double sum2 = 0.0;
			double sum1 = 0.0;

			for (int segindex = 0; segindex < m_mv.segcount; ++segindex)
			{
				stringstream filename_normals;
				stringstream filename_cloud;

				filename_normals << "normals/segment_" << segindex << ".pcd";
				filename_cloud << "segments/segment_" << segindex << ".pcd";

				// curvature point cloud
				pcl::PointCloud<pcl::Normal>::Ptr cloudc(new pcl::PointCloud<pcl::Normal>);
				if (pcl::io::loadPCDFile <pcl::Normal>(filename_normals.str(), *cloudc) == -1)
				{
					std::cout << "Cloud reading failed." << std::endl;
				}
				else
				{
					//rgb_point cloud
					pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloudrgb(new pcl::PointCloud<pcl::PointXYZRGB>);
					if (pcl::io::loadPCDFile <pcl::PointXYZRGB>(filename_cloud.str(), *cloudrgb) == -1)
					{
						std::cout << "Cloud reading failed." << std::endl;
					}
					else
					{
						double r_average = region_average(cloudc, cloudc->points.size());
						double rgb_average = region_averagergb(cloudrgb, cloudrgb->points.size());

						double heterogeneity = compute_heterogeneity(b_average, r_average);
						double rgb_heterogeneity = compute_heterogeneityrgb(b_averagergb, rgb_average);

						//std::cout<<"heterogeneity: "<<heterogeneity<<" "<<"rgb heterogeneity: "<<rgb_heterogeneity<<std::endl;

						sum2 += heterogeneity;
						sum1 += rgb_heterogeneity;

						/*std::ofstream outfile;

						outfile.open("test.txt", std::ios_base::app);
						outfile << heterogeneity; */
					}
				}
			}

			m_metrics.RGBHeterogeneity = sum1 / m_mv.segcount;

			cout << "RGB Heterogeneity: " << m_metrics.RGBHeterogeneity << endl;
		}

	private:
		int m_knn; // neighborhood size

		float m_smoothness_threshold; // smoothness threshold ... in degrees
		float m_curvature_threshold; // curvature threshold ... in degress

		int m_min_cluster; // minimum cluster size
		int m_max_cluster; // maximum cluster size

		SegmentationMetrics m_metrics;
		MetricsVariables m_mv;
	};

}

#endif // REGIONGROWINGSEGMENTATION_H
