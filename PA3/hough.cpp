#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <math.h>
#include <vector>
#include <opencv2/opencv.hpp>

#define PI 3.14
struct Node{
  int x;
  int y;
  Node* next;
};

class LinkedList{
private:
  Node* head;
  Node* ptr;
  int length;
public:
  LinkedList(){
    head = new Node();
    head->next = nullptr;
    ptr = head;
    length = 0;
  };
  LinkedList(int x,int y){
    head = new Node();
    head->x = x;
    head->y = y;
    head->next = nullptr;
    ptr = head;
    length = 1;
  }
  void addNode(int x, int y){
    if(length == 0){
      head->x = x;
      head->y = y;
      length++;
    }
    else{
    ptr->next = new Node();
    ptr = ptr->next;
    ptr->x = x;
    ptr->y = y;
    ptr->next = nullptr;
    length++;
    }
  }
  Node* getNode(){
    return ptr;
  }

  Node* getHead(){
    return head;
  }
  int getLength(){
    return length;
  }

  void clear(){
    length = 0;
  }
};

int convert(int omax, int omin, int nmax, int nmin, int value){
  int orange = (omax - omin);
  int nrange = (nmax - nmin);
  return (((value- omin) * nrange) / orange) + nmin;
}

int* hough(cv::Mat img, int threshold, int r_min, int r_max, int a_max, int b_max, int bin){
  //std::cout<< "here" << std::endl;

  int a, b;
  int new_a, new_b, new_r;
  int x_size = img.rows;
  int y_size = img.cols;
  LinkedList ll[x_size/bin][y_size/bin][r_max/bin];
  std::cout << a_max << " " << a_max/bin << std::endl;
  for(int i = 0; i < x_size; i++){
    for(int j = 0; j < y_size; j++){
      cv::Scalar val = img.at<uchar>(i,j);
      int v = abs(val[0]);
      if(v > threshold){
        std::cout << i << " " << j << std::endl;
      for(int r = r_min; r < r_max; r++){
        for(int theta = 0; theta < 360; theta++){


          a = i - r*cos(theta*PI/180);
          b = j - r*sin(theta*PI/180);
          if(a>=0 && b >=0 && a<x_size && b<y_size){
              new_a = convert(x_size,0,x_size/bin,0,a);
              new_b = convert(y_size,0,y_size/bin,0,b);
              new_r = convert(r_max,0,r_max/bin,0,r);
              //std::cout << a << " " << new_a << std::endl;
              ll[new_a][new_b][new_r].addNode(i,j);
              //std::cout << i << " " << j << std::endl;
          }
        }
        }
      }
    }
  }
  int tx, ty, tr;
  int length_max = 0;
  int* coord = new int[3];
  std::cout << x_size << " " << y_size << " " << r_max  << " " << r_min << std::endl;
  for(int i = 0; i < x_size/bin; i++){
    for(int j = 0; j < y_size/bin; j++){
      for(int r = r_min/bin; r  <r_max/bin; r++){
        if(length_max < ll[i][j][r].getLength()){

          length_max = ll[i][j][r].getLength();
          coord[0] =convert(x_size/bin,0,x_size,0,i);
          coord[1] =convert(y_size/bin,0,y_size,0,j);
          coord[2] =convert(r_max/bin,0,r_max,0,r);
          std::cout << length_max << " " << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;
        }
      }
    }
  }
  return coord;
}

int main(){
  cv::Mat src;
  cv::Mat color_src;
  cv::Mat grad_x, grad_y, abs_grad_x, abs_grad_y, grad;
  char name[20] = "overlap2.pgm";
  char* window_name = "Hough";
  src = cv::imread(name,cv::IMREAD_GRAYSCALE);
  color_src = cv::imread(name);
  cv::GaussianBlur(src,src,cv::Size(3,3),0,0,cv::BORDER_DEFAULT);
  cv::namedWindow( window_name, cv::WINDOW_AUTOSIZE );
  cv::Sobel(src,grad_x,CV_16S,1,0,3,1,0,cv::BORDER_DEFAULT);
  cv::Sobel(src, grad_y,CV_16S,0,1,3,1,0,cv::BORDER_DEFAULT);
  cv::convertScaleAbs(grad_x,abs_grad_x);
  cv::convertScaleAbs(grad_y,abs_grad_y);
  cv::addWeighted( abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad );
  //cv::imshow( window_name, grad );
  int x_size = grad.rows;
  int y_size = grad.cols;
  int length_t = 2000;
  int threshold = 50;
  int r_max = std::max(x_size,y_size);
  r_max /= 2;
  //coord = hough(grad,150, 30,100, rows, cols, 4);
  int bin = 1;
  int r_min = 30;
  int a, b;
  int new_a, new_b, new_r;

  static LinkedList ll[160][120][80];
  for(int i = 0; i < x_size; i++){
    for(int j = 0; j < y_size; j++){
      cv::Scalar val = grad.at<uchar>(i,j);
      int v = abs(val[0]);
      if(v >= threshold){
      //  std::cout << i << " " << j << std::endl;
      for(int r = r_min; r < r_max; r++){
        for(int theta = 0; theta < 360; theta++){


          a = i - r*cos(theta*PI/180);
          b = j - r*sin(theta*PI/180);
          if(a>=0 && b >=0 && a<x_size && b<y_size){
              new_a = convert(x_size,0,x_size/bin,0,a);
              new_b = convert(y_size,0,y_size/bin,0,b);
              new_r = convert(r_max,0,r_max/bin,0,r);
              //std::cout << a << " " << new_a << std::endl;
              ll[new_a][new_b][new_r].addNode(i,j);
              //std::cout << i << " " << j << std::endl;
          }
        }
        }
      }
    }
  }
  int tx, ty, tr, l = 0;
  //int* coord = new int[4];
  //std::vector<int*>keep;
  int keep [8000][4];
  for(int i = 0; i < 50; i++){
    for(int j = 0; j < 4; j++){
      keep[i][j] = 0;
    }

  }
  for(int i = 0; i < x_size/bin; i++){
    for(int j = 0; j < y_size/bin; j++){
      for(int r = r_min/bin; r  <r_max/bin; r++){
        if(length_t < ll[i][j][r].getLength()){
          //std::cout << ll[i][j][r].getLength() << " " << i << " " << j << std::endl;
          keep[l][0] = i;//convert(x_size/bin,0,x_size,0,i);
          keep[l][1] = j;//convert(y_size/bin,0,y_size,0,j);
          keep[l][2] = r;//convert(r_max/bin,0,r_max,0,r);
          keep[l][3] = ll[i][j][r].getLength();
        //  keep.push_back(coord);
          //std::cout << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;
          std::cout << keep[l][0] << " " << keep[l][1] << " " << keep[l][2] << std::endl;
          l++;
          //cv::Point center(ty,tx);
          //cv::circle(grad,center,tr,(255));
        }
        else{
          ll[i][j][r].clear();
        }
      }
    }
  }
  std::cout << std::endl;
  tx = keep[0][0];
  ty = keep[0][1];
  tr = keep[0][2];
  l = keep[0][3];
  for(int t = 1; t < 50; t++){
  //  std::cout << keep[t][0] << " " << keep[t][1] << " " << keep[t][2] << std::endl;
    //std::cout << tx << " " << ty << " " << tr << std::endl << std::endl;
    if(tx == keep[t][0] || tx - 1 == keep[t][0] || tx + 1 == keep[t][0] || ty == keep[t][1] || ty - 1 == keep[t][1] || ty + 1 == keep[t][1]){
      //std::cout << keep[t][0] << " " << keep[t][1] << " " << keep[t][2] << std::endl;
      if(l >= keep[t][3]){
        std::cout << keep[t-1][3] << " " << keep[t][3] << std::endl;
        keep[t][3] = 0;
      }
      else{
        keep[t-1][3] = 0;
        tx = keep[t][0];
        ty = keep[t][1];
        tr = keep[t][2];
        l = keep[t][3];
      }
    }
    else{
      tx = keep[t][0];
      ty = keep[t][1];
      tr = keep[t][2];
      l = keep[t][3];
    }
  }
  for(int t = 0; t < 50; t++){
    double e = 0;
    if(length_t < keep[t][3]){
      std::cout << keep[t][0] << " " << keep[t][1] << " " << keep[t][2]<< " " << keep[t][3] << std::endl;
      tx = convert(x_size/bin,0,x_size,0,keep[t][0]);
      ty = convert(y_size/bin,0,y_size,0,keep[t][1]);
      tr = convert(r_max/bin,0,r_max,0,keep[t][2]);
      std::cout << tx << " " << ty << " " << tr << std::endl;
      Node* track = ll[keep[t][0]][keep[t][1]][keep[t][2]].getHead();
      for(int i = 0; i < keep[t][3]; i++){
        int x_err = track->x;
        int y_err = track->y;
        double di = abs(sqrt(((x_err-tx)*(x_err-tx))+((y_err-ty)*(y_err-ty)))-tr);
        //std::cout << di << std::endl;
        e += di;
        track = track->next;
      }
      std::cout << e << std::endl;
      e /= keep[t][3];
      std::cout << e << std::endl;
      if(e <= .48){
        cv::Point center(ty,tx);
        cv::circle(color_src,center,tr,(255));
      }
    }
  }
  //std::cout << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;
//  cv::Point center(coord[1],coord[0]);
//  cv::circle(grad,center,coord[2],(255));
  //cv::imwrite("overlap2_circle.jpg",color_src);
  cv::imshow(window_name, color_src);
  cv::waitKey(0);
  return 0;
}
