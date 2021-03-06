// Copyright 2013, Ji Zhang, Carnegie Mellon University
// Further contributions copyright (c) 2016, Southwest Research Institute
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from this
//    software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// This is an implementation of the algorithm described in the following paper:
//   J. Zhang and S. Singh. LOAM: Lidar Odometry and Mapping in Real-time.
//     Robotics: Science and Systems Conference (RSS). Berkeley, CA, July 2014.

#include <cmath>
#include <vector>

#include <nav_msgs/Odometry.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/PointCloud2.h>
#include <tf/transform_datatypes.h>
#include <tf/transform_broadcaster.h>
#include "loam_velodyne/common.h"
#include "math_utils.h"
#include "std_msgs/Int32.h"

using std::sin;
using std::cos;
using std::atan2;

// #define M_PI   3.141592653 
#define DEG_TO_RAD(a) (a/180.0*M_PI) 
#define RAD_TO_DEG(a) (a/M_PI*180.0)

int flag_restart_scanRegistration = 0;

//added by lichunjing 2017-12-25
ros::Publisher pubRcvFlagRestartConfirmed;

double time_laser = 0.0;    

/** Point label options. */
enum PointLabel {
    CORNER_SHARP = 2,       ///< sharp corner point
    CORNER_LESS_SHARP = 1,  ///< less sharp corner point
    SURFACE_LESS_FLAT = 0,  ///< less flat surface point
    SURFACE_FLAT = -1       ///< flat surface point
};


const double scanPeriod = 0.1;

/** The number of (equally sized) regions used to distribute the feature extraction within a scan. */
int N_FEATURE_REGIONS = 6;

/** The number of surrounding points (+/- region around a point) used to calculate a point curvature. */
int CURVATURE_REGION = 5;

/** The maximum number of sharp corner points per feature region. */
int MAX_CORNER_SHARP = 2;

/** The maximum number of less sharp corner points per feature region. */
int MAX_CORNER_LESS_SHARP;

/** The maximum number of flat surface points per feature region. */
int MAX_SURFACE_FLAT = 4;

/** The voxel size used for down sizing the remaining less flat surface points. */
float LESS_FLAT_FILTER_SIZE = 0.2;

/** The curvature threshold below / above which a point is considered a flat / corner point. */
float SURFACE_CURVATURE_THRESHOLD = 0.1;

const int systemDelay = 20;
int systemInitCount = 0;
bool systemInited = false;

const int N_SCANS = 16;

float cloudCurvature[40000];
int cloudSortInd[40000];
int cloudNeighborPicked[40000];
PointLabel cloudLabel[40000];

int imuPointerFront = 0;
int imuPointerLast = -1;
const int imuQueLength = 200;

Angle imuRollStart, imuPitchStart, imuYawStart;
Angle imuRollCur, imuPitchCur, imuYawCur;

Vector3 imuVeloStart;
Vector3 imuShiftStart;
Vector3 imuVeloCur;
Vector3 imuShiftCur;

Vector3 imuShiftFromStartCur;
Vector3 imuVeloFromStartCur;

double imuTime[imuQueLength] = {0};
float imuRoll[imuQueLength] = {0};
float imuPitch[imuQueLength] = {0};
float imuYaw[imuQueLength] = {0};

Vector3 imuAcc[imuQueLength];
Vector3 imuVelo[imuQueLength];
Vector3 imuShift[imuQueLength];

ros::Publisher pubLaserCloud;
ros::Publisher pubCornerPointsSharp;
ros::Publisher pubCornerPointsLessSharp;
ros::Publisher pubSurfPointsFlat;
ros::Publisher pubSurfPointsLessFlat;
ros::Publisher pubImuTrans;

void ShiftToStartIMU(float pointTime)
{
  imuShiftFromStartCur = imuShiftCur - imuShiftStart - imuVeloStart * pointTime;

  Vector3 v1 = rotateY( imuShiftFromStartCur, -imuYawStart);
  Vector3 v2 = rotateX( v1,                   -imuPitchStart);
  imuShiftFromStartCur = rotateZ( v2,         -imuRollStart);
}

void VeloToStartIMU()
{
  imuVeloFromStartCur = imuVeloCur - imuVeloStart;

  Vector3 v1 = rotateY( imuVeloFromStartCur, -imuYawStart);
  Vector3 v2 = rotateX( v1, -imuPitchStart);
  imuVeloFromStartCur = rotateZ( v2, -imuRollStart);
}

void TransformToStartIMU(PointType *p)
{
  Vector3 v1 = rotateZ( *p, imuRollCur);
  Vector3 v2 = rotateX( v1, imuPitchCur);
  Vector3 v3 = rotateY( v2, imuYawCur);

  Vector3 v4 = rotateY( v3, -imuYawStart);
  Vector3 v5 = rotateX( v4, -imuPitchStart);
  Vector3 v6 = rotateZ( v5, -imuRollStart);

  v6 += imuShiftFromStartCur;

  p->x = v6.x();
  p->y = v6.y();
  p->z = v6.z();
}

void AccumulateIMUShift()
{
  float roll = imuRoll[imuPointerLast];
  float pitch = imuPitch[imuPointerLast];
  float yaw = imuYaw[imuPointerLast];
  Vector3 acc = imuAcc[imuPointerLast];

  Vector3 v1 = rotateZ( acc, roll );
  Vector3 v2 = rotateX( v1, pitch );
  acc        = rotateY( v2, yaw );


  int imuPointerBack = (imuPointerLast + imuQueLength - 1) % imuQueLength;
  double timeDiff = imuTime[imuPointerLast] - imuTime[imuPointerBack];
  if (timeDiff < scanPeriod) {

    imuShift[imuPointerLast] = imuShift[imuPointerBack] + (imuVelo[imuPointerBack] * timeDiff)
                              + acc * (0.5* timeDiff * timeDiff);

    imuVelo[imuPointerLast] = imuVelo[imuPointerBack] + acc * timeDiff;
  }
}

void laserCloudHandler(const sensor_msgs::PointCloud2ConstPtr& laserCloudMsg)
{
  if (!systemInited) {
    systemInitCount++;
    if (systemInitCount >= systemDelay) {
      systemInited = true;
    }
    return;
  }

  std::vector<int> scanStartInd(N_SCANS, 0);
  std::vector<int> scanEndInd(N_SCANS, 0);
  
  double timeScanCur = laserCloudMsg->header.stamp.toSec();
  pcl::PointCloud<pcl::PointXYZ> laserCloudIn;
  pcl::fromROSMsg(*laserCloudMsg, laserCloudIn);
  std::vector<int> indices;
  pcl::removeNaNFromPointCloud(laserCloudIn, laserCloudIn, indices);
  int cloudSize = laserCloudIn.points.size();
  float startOri = -atan2(laserCloudIn.points[0].y, laserCloudIn.points[0].x);
  float endOri = -atan2(laserCloudIn.points[cloudSize - 1].y,
                        laserCloudIn.points[cloudSize - 1].x) + 2 * M_PI;

  time_laser = timeScanCur;
  // ROS_INFO("laser time:%f",time_laser);

  if (endOri - startOri > 3 * M_PI) {
    endOri -= 2 * M_PI;
  } else if (endOri - startOri < M_PI) {
    endOri += 2 * M_PI;
  }
  bool halfPassed = false;
  int count = cloudSize;
  PointType point;
  std::vector<pcl::PointCloud<PointType> > laserCloudScans(N_SCANS);
  for (int i = 0; i < cloudSize; i++) {
    point.x = laserCloudIn.points[i].y;
    point.y = laserCloudIn.points[i].z;
    point.z = laserCloudIn.points[i].x;

    float angle = atan(point.y / sqrt(point.x * point.x + point.z * point.z)) * 180 / M_PI;
    int scanID;
    int roundedAngle = int(angle + (angle<0.0?-0.5:+0.5)); 
    if (roundedAngle > 0){
      scanID = roundedAngle;
    }
    else {
      scanID = roundedAngle + (N_SCANS - 1);
    }
    if (scanID > (N_SCANS - 1) || scanID < 0 ){
      count--;
      continue;
    }

    float ori = -atan2(point.x, point.z);
    if (!halfPassed) {
      if (ori < startOri - M_PI / 2) {
        ori += 2 * M_PI;
      } else if (ori > startOri + M_PI * 3 / 2) {
        ori -= 2 * M_PI;
      }

      if (ori - startOri > M_PI) {
        halfPassed = true;
      }
    } else {
      ori += 2 * M_PI;

      if (ori < endOri - M_PI * 3 / 2) {
        ori += 2 * M_PI;
      } else if (ori > endOri + M_PI / 2) {
        ori -= 2 * M_PI;
      } 
    }

    float relTime = (ori - startOri) / (endOri - startOri);
    point.intensity = scanID + scanPeriod * relTime;

    if (imuPointerLast >= 0) {
      float pointTime = relTime * scanPeriod;
      while (imuPointerFront != imuPointerLast) {
        if (timeScanCur + pointTime < imuTime[imuPointerFront]) {
          break;
        }
        imuPointerFront = (imuPointerFront + 1) % imuQueLength;
      }

      if (timeScanCur + pointTime > imuTime[imuPointerFront]) {
        imuRollCur = imuRoll[imuPointerFront];
        imuPitchCur = imuPitch[imuPointerFront];
        imuYawCur = imuYaw[imuPointerFront];

        imuVeloCur = imuVelo[imuPointerFront];
        imuShiftCur = imuShift[imuPointerFront];
      } else {
        int imuPointerBack = (imuPointerFront + imuQueLength - 1) % imuQueLength;
        float ratioFront = (timeScanCur + pointTime - imuTime[imuPointerBack]) 
                         / (imuTime[imuPointerFront] - imuTime[imuPointerBack]);
        float ratioBack = (imuTime[imuPointerFront] - timeScanCur - pointTime) 
                        / (imuTime[imuPointerFront] - imuTime[imuPointerBack]);

        imuRollCur = imuRoll[imuPointerFront] * ratioFront + imuRoll[imuPointerBack] * ratioBack;
        imuPitchCur = imuPitch[imuPointerFront] * ratioFront + imuPitch[imuPointerBack] * ratioBack;
        if (imuYaw[imuPointerFront] - imuYaw[imuPointerBack] > M_PI) {
          imuYawCur = imuYaw[imuPointerFront] * ratioFront + (imuYaw[imuPointerBack] + 2 * M_PI) * ratioBack;
        } else if (imuYaw[imuPointerFront] - imuYaw[imuPointerBack] < -M_PI) {
          imuYawCur = imuYaw[imuPointerFront] * ratioFront + (imuYaw[imuPointerBack] - 2 * M_PI) * ratioBack;
        } else {
          imuYawCur = imuYaw[imuPointerFront] * ratioFront + imuYaw[imuPointerBack] * ratioBack;
        }

        imuVeloCur = imuVelo[imuPointerFront] * ratioFront + imuVelo[imuPointerBack] * ratioBack;
        imuShiftCur = imuShift[imuPointerFront] * ratioFront + imuShift[imuPointerBack] * ratioBack;
      }
      if (i == 0) {
        imuRollStart = imuRollCur;
        imuPitchStart = imuPitchCur;
        imuYawStart = imuYawCur;

        imuVeloStart = imuVeloCur;
        imuShiftStart = imuShiftCur;

      } else {
        ShiftToStartIMU(pointTime);
        VeloToStartIMU();
        TransformToStartIMU(&point);
      }
    }
    laserCloudScans[scanID].push_back(point);
  }
  cloudSize = count;

  pcl::PointCloud<PointType>::Ptr laserCloud(new pcl::PointCloud<PointType>());
  for (int i = 0; i < N_SCANS; i++) {
    *laserCloud += laserCloudScans[i];
  }

  int scanCount = -1;
  float pointWeight = -2 * CURVATURE_REGION;
  for (int i = CURVATURE_REGION; i < cloudSize - CURVATURE_REGION; i++) {
    float diffX = pointWeight * laserCloud->points[i].x;
    float diffY = pointWeight * laserCloud->points[i].y;
    float diffZ = pointWeight * laserCloud->points[i].z;

    for (int j = 1; j <= CURVATURE_REGION; j++) {
      diffX += laserCloud->points[i + j].x + laserCloud->points[i - j].x;
      diffY += laserCloud->points[i + j].y + laserCloud->points[i - j].y;
      diffZ += laserCloud->points[i + j].z + laserCloud->points[i - j].z;
    }

    cloudCurvature[i] = diffX * diffX + diffY * diffY + diffZ * diffZ;
    cloudSortInd[i] = i;
    cloudNeighborPicked[i] = 0;
    cloudLabel[i] = SURFACE_LESS_FLAT;

    if (int(laserCloud->points[i].intensity) != scanCount) {
      scanCount = int(laserCloud->points[i].intensity);

      if (scanCount > 0 && scanCount < N_SCANS) {
        scanStartInd[scanCount] = i + CURVATURE_REGION;
        scanEndInd[scanCount - 1] = i - CURVATURE_REGION;
      }
    }
  }
  scanStartInd[0] = CURVATURE_REGION;
  scanEndInd.back() = cloudSize - CURVATURE_REGION;

  for (int i = CURVATURE_REGION; i < cloudSize - CURVATURE_REGION - 1; i++) {
    float diffX = laserCloud->points[i + 1].x - laserCloud->points[i].x;
    float diffY = laserCloud->points[i + 1].y - laserCloud->points[i].y;
    float diffZ = laserCloud->points[i + 1].z - laserCloud->points[i].z;
    float diff = diffX * diffX + diffY * diffY + diffZ * diffZ;

    if (diff > 0.1) {

      float depth1 = sqrt(laserCloud->points[i].x * laserCloud->points[i].x + 
                     laserCloud->points[i].y * laserCloud->points[i].y +
                     laserCloud->points[i].z * laserCloud->points[i].z);

      float depth2 = sqrt(laserCloud->points[i + 1].x * laserCloud->points[i + 1].x + 
                     laserCloud->points[i + 1].y * laserCloud->points[i + 1].y +
                     laserCloud->points[i + 1].z * laserCloud->points[i + 1].z);

      if (depth1 > depth2) {
        diffX = laserCloud->points[i + 1].x - laserCloud->points[i].x * depth2 / depth1;
        diffY = laserCloud->points[i + 1].y - laserCloud->points[i].y * depth2 / depth1;
        diffZ = laserCloud->points[i + 1].z - laserCloud->points[i].z * depth2 / depth1;

        if (sqrt(diffX * diffX + diffY * diffY + diffZ * diffZ) / depth2 < 0.1) {
          for (int j = 0; j <= CURVATURE_REGION; j++) {
            cloudNeighborPicked[i - j] = 1;
          }
        }
      } else {
        diffX = laserCloud->points[i + 1].x * depth1 / depth2 - laserCloud->points[i].x;
        diffY = laserCloud->points[i + 1].y * depth1 / depth2 - laserCloud->points[i].y;
        diffZ = laserCloud->points[i + 1].z * depth1 / depth2 - laserCloud->points[i].z;

        if (sqrt(diffX * diffX + diffY * diffY + diffZ * diffZ) / depth1 < 0.1) {
          for (int j = CURVATURE_REGION + 1; j > 0 ; j--) {
            cloudNeighborPicked[i + j] = 1;
          }
        }
      }
    }

    float diffX2 = laserCloud->points[i].x - laserCloud->points[i - 1].x;
    float diffY2 = laserCloud->points[i].y - laserCloud->points[i - 1].y;
    float diffZ2 = laserCloud->points[i].z - laserCloud->points[i - 1].z;
    float diff2 = diffX2 * diffX2 + diffY2 * diffY2 + diffZ2 * diffZ2;

    float dis = laserCloud->points[i].x * laserCloud->points[i].x
              + laserCloud->points[i].y * laserCloud->points[i].y
              + laserCloud->points[i].z * laserCloud->points[i].z;

    if (diff > 0.0002 * dis && diff2 > 0.0002 * dis) {
      cloudNeighborPicked[i] = 1;
    }
  }


  pcl::PointCloud<PointType> cornerPointsSharp;
  pcl::PointCloud<PointType> cornerPointsLessSharp;
  pcl::PointCloud<PointType> surfPointsFlat;
  pcl::PointCloud<PointType> surfPointsLessFlat;

  for (int i = 0; i < N_SCANS; i++) {
    pcl::PointCloud<PointType>::Ptr surfPointsLessFlatScan(new pcl::PointCloud<PointType>);
    for (int j = 0; j < N_FEATURE_REGIONS; j++) {
      int sp = (scanStartInd[i] * (N_FEATURE_REGIONS - j)  + scanEndInd[i] * j) / N_FEATURE_REGIONS;
      int ep = (scanStartInd[i] * (N_FEATURE_REGIONS - 1 - j)  + scanEndInd[i] * (j + 1)) / N_FEATURE_REGIONS - 1;

      for (int k = sp + 1; k <= ep; k++) {
        for (int l = k; l >= sp + 1; l--) {
          if (cloudCurvature[cloudSortInd[l]] < cloudCurvature[cloudSortInd[l - 1]]) {
            int temp = cloudSortInd[l - 1];
            cloudSortInd[l - 1] = cloudSortInd[l];
            cloudSortInd[l] = temp;
          }
        }
      }

      int largestPickedNum = 0;
      for (int k = ep; k >= sp; k--) {
        int ind = cloudSortInd[k];
        if (cloudNeighborPicked[ind] == 0 &&
            cloudCurvature[ind] > SURFACE_CURVATURE_THRESHOLD) {

          largestPickedNum++;
          if (largestPickedNum <= MAX_CORNER_SHARP) {
            cloudLabel[ind] = CORNER_SHARP;
            cornerPointsSharp.push_back(laserCloud->points[ind]);
            cornerPointsLessSharp.push_back(laserCloud->points[ind]);
          } else if (largestPickedNum <= MAX_CORNER_LESS_SHARP) {
            cloudLabel[ind] = CORNER_LESS_SHARP;
            cornerPointsLessSharp.push_back(laserCloud->points[ind]);
          } else {
            break;
          }

          cloudNeighborPicked[ind] = 1;
          for (int l = 1; l <= CURVATURE_REGION; l++) {
            float diffX = laserCloud->points[ind + l].x 
                        - laserCloud->points[ind + l - 1].x;
            float diffY = laserCloud->points[ind + l].y 
                        - laserCloud->points[ind + l - 1].y;
            float diffZ = laserCloud->points[ind + l].z 
                        - laserCloud->points[ind + l - 1].z;
            if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05) {
              break;
            }

            cloudNeighborPicked[ind + l] = 1;
          }
          for (int l = -1; l >= -CURVATURE_REGION; l--) {
            float diffX = laserCloud->points[ind + l].x 
                        - laserCloud->points[ind + l + 1].x;
            float diffY = laserCloud->points[ind + l].y 
                        - laserCloud->points[ind + l + 1].y;
            float diffZ = laserCloud->points[ind + l].z 
                        - laserCloud->points[ind + l + 1].z;
            if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05) {
              break;
            }

            cloudNeighborPicked[ind + l] = 1;
          }
        }
      }

      int smallestPickedNum = 0;
      for (int k = sp; k <= ep; k++) {
        int ind = cloudSortInd[k];
        if (cloudNeighborPicked[ind] == 0 &&
            cloudCurvature[ind] < SURFACE_CURVATURE_THRESHOLD) {

          cloudLabel[ind] = SURFACE_FLAT;
          surfPointsFlat.push_back(laserCloud->points[ind]);

          smallestPickedNum++;
          if (smallestPickedNum >= MAX_SURFACE_FLAT) {
            break;
          }

          cloudNeighborPicked[ind] = 1;
          for (int l = 1; l <= CURVATURE_REGION; l++) {
            float diffX = laserCloud->points[ind + l].x 
                        - laserCloud->points[ind + l - 1].x;
            float diffY = laserCloud->points[ind + l].y 
                        - laserCloud->points[ind + l - 1].y;
            float diffZ = laserCloud->points[ind + l].z 
                        - laserCloud->points[ind + l - 1].z;
            if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05) {
              break;
            }

            cloudNeighborPicked[ind + l] = 1;
          }
          for (int l = -1; l >= -CURVATURE_REGION; l--) {
            float diffX = laserCloud->points[ind + l].x 
                        - laserCloud->points[ind + l + 1].x;
            float diffY = laserCloud->points[ind + l].y 
                        - laserCloud->points[ind + l + 1].y;
            float diffZ = laserCloud->points[ind + l].z 
                        - laserCloud->points[ind + l + 1].z;
            if (diffX * diffX + diffY * diffY + diffZ * diffZ > 0.05) {
              break;
            }

            cloudNeighborPicked[ind + l] = 1;
          }
        }
      }

      for (int k = sp; k <= ep; k++) {
        if (cloudLabel[k] <= SURFACE_LESS_FLAT) {
          surfPointsLessFlatScan->push_back(laserCloud->points[k]);
        }
      }
    }

    pcl::PointCloud<PointType> surfPointsLessFlatScanDS;
    pcl::VoxelGrid<PointType> downSizeFilter;
    downSizeFilter.setInputCloud(surfPointsLessFlatScan);
    downSizeFilter.setLeafSize(LESS_FLAT_FILTER_SIZE, LESS_FLAT_FILTER_SIZE, LESS_FLAT_FILTER_SIZE);
    downSizeFilter.filter(surfPointsLessFlatScanDS);

    surfPointsLessFlat += surfPointsLessFlatScanDS;
  }

  sensor_msgs::PointCloud2 laserCloudOutMsg;
  pcl::toROSMsg(*laserCloud, laserCloudOutMsg);
  laserCloudOutMsg.header.stamp = laserCloudMsg->header.stamp;
  laserCloudOutMsg.header.frame_id = "/camera";
  pubLaserCloud.publish(laserCloudOutMsg);

  sensor_msgs::PointCloud2 cornerPointsSharpMsg;
  pcl::toROSMsg(cornerPointsSharp, cornerPointsSharpMsg);
  cornerPointsSharpMsg.header.stamp = laserCloudMsg->header.stamp;
  cornerPointsSharpMsg.header.frame_id = "/camera";
  pubCornerPointsSharp.publish(cornerPointsSharpMsg);

  sensor_msgs::PointCloud2 cornerPointsLessSharpMsg;
  pcl::toROSMsg(cornerPointsLessSharp, cornerPointsLessSharpMsg);
  cornerPointsLessSharpMsg.header.stamp = laserCloudMsg->header.stamp;
  cornerPointsLessSharpMsg.header.frame_id = "/camera";
  pubCornerPointsLessSharp.publish(cornerPointsLessSharpMsg);

  sensor_msgs::PointCloud2 surfPointsFlat2;
  pcl::toROSMsg(surfPointsFlat, surfPointsFlat2);
  surfPointsFlat2.header.stamp = laserCloudMsg->header.stamp;
  surfPointsFlat2.header.frame_id = "/camera";
  pubSurfPointsFlat.publish(surfPointsFlat2);

  sensor_msgs::PointCloud2 surfPointsLessFlat2;
  pcl::toROSMsg(surfPointsLessFlat, surfPointsLessFlat2);
  surfPointsLessFlat2.header.stamp = laserCloudMsg->header.stamp;
  surfPointsLessFlat2.header.frame_id = "/camera";
  pubSurfPointsLessFlat.publish(surfPointsLessFlat2);

  pcl::PointCloud<pcl::PointXYZ> imuTrans(4, 1);
  imuTrans.points[0].x = imuPitchStart.value();
  imuTrans.points[0].y = imuYawStart.value();
  imuTrans.points[0].z = imuRollStart.value();

  imuTrans.points[1].x = imuPitchCur.value();
  imuTrans.points[1].y = imuYawCur.value();
  imuTrans.points[1].z = imuRollCur.value();

  imuTrans.points[2].x = imuShiftFromStartCur.x();
  imuTrans.points[2].y = imuShiftFromStartCur.y();
  imuTrans.points[2].z = imuShiftFromStartCur.z();

  imuTrans.points[3].x = imuVeloFromStartCur.x();
  imuTrans.points[3].y = imuVeloFromStartCur.y();
  imuTrans.points[3].z = imuVeloFromStartCur.z();

  sensor_msgs::PointCloud2 imuTransMsg;
  pcl::toROSMsg(imuTrans, imuTransMsg);
  imuTransMsg.header.stamp = laserCloudMsg->header.stamp;
  imuTransMsg.header.frame_id = "/camera";
  pubImuTrans.publish(imuTransMsg);
}

void imuHandler(const sensor_msgs::Imu::ConstPtr& imuIn)
{
  double roll, pitch, yaw;
  tf::Quaternion orientation;
  tf::quaternionMsgToTF(imuIn->orientation, orientation);
  tf::Matrix3x3(orientation).getRPY(roll, pitch, yaw);
  //tf::Matrix3x3(orientation).getRPY(pitch, yaw, roll);

  Vector3 acc;
  acc.x() = imuIn->linear_acceleration.y - sin(roll) * cos(pitch) * 9.81;
  acc.y() = imuIn->linear_acceleration.z - cos(roll) * cos(pitch) * 9.81;
  acc.z() = imuIn->linear_acceleration.x + sin(pitch) * 9.81;

  // acc.x() = 0.0;
  // acc.y() = 0.0;
  // acc.z() = 0.0;

  imuPointerLast = (imuPointerLast + 1) % imuQueLength;

  imuTime[imuPointerLast] = imuIn->header.stamp.toSec();
  imuRoll[imuPointerLast] = roll;
  imuPitch[imuPointerLast] = pitch;
  imuYaw[imuPointerLast] = yaw;
  imuAcc[imuPointerLast] = acc;

  // ROS_INFO("acc_x=%f  acc_y=%f  acc_z=%f  ",acc.x(), acc.y(), acc.z());
  // ROS_INFO("Quaternion:%.3f %.3f %.3f %.3f", imuIn->orientation.x, imuIn->orientation.y, imuIn->orientation.z, imuIn->orientation.w);
  // ROS_INFO("laser time:%f  Acc:%f %f %f  Angle:%f %f %f",time_laser, imuIn->linear_acceleration.x, imuIn->linear_acceleration.y, imuIn->linear_acceleration.z, RAD_TO_DEG(roll), RAD_TO_DEG(pitch), RAD_TO_DEG(yaw));
  // ROS_INFO("laser time:%f  imu time:%f",time_laser, imuIn->header.stamp.toSec());

  AccumulateIMUShift();
}

//added by lichunjing 2017-12-26
void RestartscanRegistrationHandler(const std_msgs::Int32::ConstPtr& msgIn)
{
  std_msgs::Int32 msg;
  msg.data = 1;

  if(msgIn->data == 1)
  {
    flag_restart_scanRegistration = 1;
    pubRcvFlagRestartConfirmed.publish(msg);
  }

  if((msgIn->data == 0)&&(flag_restart_scanRegistration == 1))
  {
    ROS_FATAL("restart_loam: scanRegistration");
    ros::shutdown();
  }
}

int main(int argc, char** argv)
{
  ros::init(argc, argv, "scanRegistration");
  ros::NodeHandle nh;


  N_FEATURE_REGIONS = nh.param("/scanRegistration/featureRegions", N_FEATURE_REGIONS);
  if (N_FEATURE_REGIONS < 1) {
    ROS_FATAL("Invalid featureRegions parameter: %d (expected >= 1)", N_FEATURE_REGIONS);
    ros::shutdown();
  }

  CURVATURE_REGION = nh.param("/scanRegistration/curvatureRegion", CURVATURE_REGION);
  if (CURVATURE_REGION < 1) {
    ROS_FATAL("Invalid curvatureRegion parameter: %d (expected >= 1)", CURVATURE_REGION);
    ros::shutdown();
  }

  MAX_CORNER_SHARP = nh.param("/scanRegistration/maxCornerSharp", MAX_CORNER_SHARP);
  if (MAX_CORNER_SHARP < 1) {
    ROS_FATAL("Invalid maxCornerSharp parameter: %d (expected >= 1)", MAX_CORNER_SHARP);
    ros::shutdown();
  }

  MAX_CORNER_LESS_SHARP = nh.param("/scanRegistration/maxCornerLessSharp", 10 * MAX_CORNER_SHARP);
  if (MAX_CORNER_LESS_SHARP < MAX_CORNER_SHARP) {
    ROS_FATAL("Invalid maxCornerLessSharp parameter: %d (expected >= %d)", MAX_CORNER_LESS_SHARP, MAX_CORNER_SHARP);
    ros::shutdown();
  }

  MAX_SURFACE_FLAT = nh.param("/scanRegistration/maxSurfaceFlat", MAX_SURFACE_FLAT);
  if (MAX_SURFACE_FLAT < 1) {
    ROS_FATAL("Invalid maxSurfaceFlat parameter: %d (expected >= 1)", MAX_SURFACE_FLAT);
    ros::shutdown();
  }

  SURFACE_CURVATURE_THRESHOLD = nh.param("/scanRegistration/surfaceCurvatureThreshold", SURFACE_CURVATURE_THRESHOLD);
  if (SURFACE_CURVATURE_THRESHOLD < 0.001) {
    ROS_FATAL("Invalid surfaceCurvatureThreshold parameter: %f (expected >= 0.001)", SURFACE_CURVATURE_THRESHOLD);
    ros::shutdown();
  }

  LESS_FLAT_FILTER_SIZE = nh.param("/scanRegistration/lessFlatFilterSize", LESS_FLAT_FILTER_SIZE);
  if (LESS_FLAT_FILTER_SIZE < 0.001) {
    ROS_FATAL("Invalid lessFlatFilterSize parameter: %f (expected >= 0.001)", LESS_FLAT_FILTER_SIZE);
    ros::shutdown();
  }


  ROS_INFO("Using  %d  feature regions per scan.", N_FEATURE_REGIONS);
  ROS_INFO("Using  +/- %d  points for curvature calculation.", CURVATURE_REGION);
  ROS_INFO("Using at most  %d  sharp  and  %d  less sharp corner points per feature region.",
           MAX_CORNER_SHARP, MAX_CORNER_LESS_SHARP);
  ROS_INFO("Using at most  %d  flat surface points per feature region.", MAX_SURFACE_FLAT);
  ROS_INFO("Using  %g  as surface curvature threshold.", SURFACE_CURVATURE_THRESHOLD);
  ROS_INFO("Using  %g  as less flat surface points voxel filter size.", LESS_FLAT_FILTER_SIZE);


  ros::Subscriber subLaserCloud = nh.subscribe<sensor_msgs::PointCloud2> 
                                  ("/velodyne_points", 2, laserCloudHandler);

  ros::Subscriber subImu = nh.subscribe<sensor_msgs::Imu> ("/imu/data", 50, imuHandler);

  //added by lichunjing 2017-12-26
  ros::Subscriber subscanRegistration = nh.subscribe<std_msgs::Int32>("/flag_restart_scanRegistration", 20, RestartscanRegistrationHandler);

  pubLaserCloud = nh.advertise<sensor_msgs::PointCloud2>
                                 ("/velodyne_cloud_2", 2);

  pubCornerPointsSharp = nh.advertise<sensor_msgs::PointCloud2>
                                        ("/laser_cloud_sharp", 2);

  pubCornerPointsLessSharp = nh.advertise<sensor_msgs::PointCloud2>
                                            ("/laser_cloud_less_sharp", 2);

  pubSurfPointsFlat = nh.advertise<sensor_msgs::PointCloud2>
                                       ("/laser_cloud_flat", 2);

  pubSurfPointsLessFlat = nh.advertise<sensor_msgs::PointCloud2>
                                           ("/laser_cloud_less_flat", 2);

  pubImuTrans = nh.advertise<sensor_msgs::PointCloud2> ("/imu_trans", 5);

  //added by lichunjing 2017-12-26
  pubRcvFlagRestartConfirmed = nh.advertise<std_msgs::Int32>("/flag_rcved_scanRegistration", 2);

  ros::spin();

  return 0;
}

