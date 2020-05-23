#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <opencv2/opencv.hpp>

double error(double r, double r0, double c, double c0){
  return sqrt(((r-r0)*(r-r0)) + ((c-c0)*(c-c0)));
}

cv::Point2f Distortion(cv::Mat d, double fx, double fy, double ox, double oy,cv::Point2f coord){
  std::cout << d << std::endl;
  double r2 = coord.x * coord.x + coord.y * coord.y;
  double r = sqrt(r2);
  double r4 = r * r * r * r;
  double r6 = r * r * r * r * r * r;
  cv::Scalar k1 = d.at<double>(0);
  cv::Scalar k2 = d.at<double>(1);
  cv::Scalar k3 = d.at<double>(4);
  cv::Scalar p1 = d.at<double>(2);
  cv::Scalar p2 = d.at<double>(3);
  double k_1 = k1[0];
  double k_2 = k2[0];
  double k_3 = k3[0];
  double p_1 = p1[0];
  double p_2 = p2[0];
  std::cout << k_1 << " " << k_2 << " " << p_1 << " " << p_2 << " " << k_3 << std::endl;
  double tan_x = coord.x * (1 + k_1*r2 + k_2*r4 + k_3*r6);
  double tan_y = coord.y * (1 + k_1*r2 + k_2*r4 + k_3*r6);
  double rad_x = 2*p_1*coord.y+ p_2*(r2+2*coord.x*coord.x);
  double rad_y = p_1*(r2 + 2*coord.y*coord.y) + 2*p_2*coord.x;
  cv::Point2f projected(fx* (tan_x+rad_x) + ox,fy * (tan_y+rad_y) + oy);
  return projected;
}

void projectionMatrix(std::vector<cv::Point3f> points, cv::Mat rot, cv::Mat tran, cv::Mat camera, cv::Mat dist, std::vector<cv::Point2f> & project){
  cv::Mat rot_mat;

  cv::Scalar f_x = camera.at<double>(0,0);
  cv::Scalar f_y = camera.at<double>(1,1);
  cv::Scalar O_x = camera.at<double>(0,2);
  cv::Scalar O_y = camera.at<double>(1,2);

  double fx = f_x[0];
  double fy = f_y[0];
  double ox = O_x[0];
  double oy = O_y[0];

  cv::Rodrigues(rot,rot_mat);

//  std::cout << rot << std::endl;
//  std::cout << rot_mat << std::endl;

  cv::Scalar r_11 = rot_mat.at<double>(0,0);
  cv::Scalar r_12 = rot_mat.at<double>(0,1);
  cv::Scalar r_13 = rot_mat.at<double>(0,2);
  cv::Scalar t_1 = tran.at<double>(0,0);

  double r11 = r_11[0];
  double r12 = r_12[0];
  double r13 = r_13[0];
  double t1 = t_1[0];

  cv::Scalar r_21 = rot_mat.at<double>(1,0);
  cv::Scalar r_22 = rot_mat.at<double>(1,1);
  cv::Scalar r_23 = rot_mat.at<double>(1,2);
  cv::Scalar t_2 = tran.at<double>(0,1);

  double r21 = r_21[0];
  double r22 = r_22[0];
  double r23 = r_23[0];
  double t2 = t_2[0];

  cv::Scalar r_31 = rot_mat.at<double>(2,0);
  cv::Scalar r_32 = rot_mat.at<double>(2,1);
  cv::Scalar r_33 = rot_mat.at<double>(2,2);
  cv::Scalar t_3 = tran.at<double>(0,2);

  double r31 = r_31[0];
  double r32 = r_32[0];
  double r33 = r_33[0];
  double t3 = t_3[0];
  //std::cout << tran << std::endl;
  //std::cout << t_1 << " " << t_2 << " " << t_3 << std::endl;
  std::cout << r11 << " " << r12 << " " << r13 << " "<<  t1 << std::endl;
  std::cout << r21 << " " << r22 << " " << r23 << " " << t2 << std::endl;
  std::cout << r31 << " " << r32 << " " << r33 << " " << t3 << std::endl;
  double projectMat[3][4] = {{fx*r11+ox*r31, fx*r12+ox*r32, fx*r13+ox*r33,fx*t1+ox*t3},
                             {fy*r21+oy*r31, fy*r22+oy*r32, fy*r23+oy*r33,fy*t2+oy*t3},
                             {r31,r32,r33,t3}};
  for(int i = 0; i < points.size(); i++){
    cv::Point3f c = points[i];
    double xim = (projectMat[0][0]*c.x + projectMat[0][1]*c.y + projectMat[0][2]*c.z + projectMat[0][3])/
                 (projectMat[2][0]*c.x + projectMat[2][1]*c.y + projectMat[2][2]*c.z + projectMat[2][3]);
    double yim = (projectMat[1][0]*c.x + projectMat[1][1]*c.y + projectMat[1][2]*c.z + projectMat[1][3])/
                 (projectMat[2][0]*c.x + projectMat[2][1]*c.y + projectMat[2][2]*c.z + projectMat[2][3]);
  //  std::cout << xim << " " << yim << std::endl;
    cv::Point2f coordinates(xim,yim);
    //coordinates = Distortion(dist,fx,fy,ox,oy,coordinates);
  //  std::cout << coordinates << std::endl;
    project.push_back(coordinates);
  }
}


int main(){
  //char[20] pixname;
  //char[20] worldname;
  std::string pixname = "./pixel/*.txt";
  std::string worldname = "./world/*.txt";
  std::string s_img = "./images/*.bmp";

  std::vector<cv::String> pixs;
  std::vector<cv::String> worlds;
  std::vector<cv::String> imgs;


  cv::glob(pixname,pixs);
  cv::glob(worldname,worlds);
  cv::glob(s_img,imgs);

  std::ifstream pixfile;
  std::ifstream worldfile;



  std::vector<std::vector<cv::Point3f>> worldpts;
  std::vector<std::vector<cv::Point2f>> pixpts;



  double px, py;
  for(int i = 0; i < 15; i++){
    std::vector<cv::Point2f> pixpt;
    std::vector<cv::Point3f> worldpt;
    pixfile.open(pixs[i]);
    std::cout << pixs[i] << std::endl;
    worldfile.open(worlds[i]);
    while(pixfile){
      pixfile >> px;
      pixfile >> py;
      //std::cout << px << " " << py << std::endl;
      pixpt.push_back(cv::Point2f(px,py));
    }
    while(worldfile){
      worldfile >> px;
      worldfile >> py;
      worldpt.push_back(cv::Point3f(px,py,0));
    }
    pixfile.close();
    worldfile.close();
    worldpts.push_back(worldpt);
    pixpts.push_back(pixpt);

  }

//  std::cout << "here" << std::endl;

  int rows = 640;
  int cols = 480;
  cv::Mat camera,dist,r,t,d;
  //std::vector<std::vector<int>> r,t;
  //cv::Mat un_img; //= cv::imread("./images/cap01.bmp");
  cv::calibrateCamera(worldpts,pixpts,cv::Size(rows,cols),camera,dist,r,t);
  std::cout << "cameraMatrix: " << camera << std::endl;
  std::cout << "dist: " << dist <<std::endl;
  std::cout << "rot: " << r << std::endl;
  std::cout << "tran: " << t << std::endl;
  double total_image_e = 0;
  for(int i = 0; i < 15; i++){
    double e = 0;
    cv::Mat img = cv::imread(imgs[i]);
    std::vector<cv::Point2f> projected;
    cv::Mat rrow = r.row(i);
    cv::Mat trow = t.row(i);
    projectionMatrix(worldpts[i],rrow,trow,camera,dist,projected);
    std::cout << "Error for image: " << i << std::endl;
    for(int j = 0; j < projected.size(); j++){
      e += error(projected[j].x,pixpts[i][j].x,projected[j].y,pixpts[i][j].y);
      cv::circle(img,projected[j],4,cv::Scalar(0,0,255),-1);
    }
      e /= projected.size();
      //std::cout << e << std::endl;
      total_image_e += e;
    std::string img_out = "image_" + std::to_string(i) + ".bmp";
    cv::imwrite(img_out,img);
  }
  total_image_e /= 15;
  //std::cout << total_image_e << std::endl;
  return 0;
}
