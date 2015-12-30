#include "opencv2/video/tracking.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/videoio/videoio.hpp"
#include "opencv2/highgui/highgui.hpp"

#include <iostream>
#include <fstream>
#include <map>
#include <list>
#include <ctype.h>
#include <thread>        
#include <mutex>          

using namespace cv;
using namespace std;

#define MAXTIME 300

Mat frame;
mutex mtxFrame;
mutex mtxProbe;
map<int,Rect> probe;
map<int,list<Vec4f> >probeValue;
Point lastPoint;
int indFrame=0;
double tpsFrame=0;
char mode = 's';
vector<Scalar> color={Scalar(255,0,0),Scalar(0,255,0),Scalar(0,0,255),Scalar(255,255,0),Scalar(255,255,255)};



int main(int ac, char** av) 
{
    Mat m = imread("C:/Users/Laurent.PC-LAURENT-VISI/Documents/Visual Studio 2013/HomomorphicFiltering/20150727_193940Homomorphic.jpg",IMREAD_COLOR );
//    Mat m = imread("f:/lib/opencv/samples/data/basketball1.png",IMREAD_COLOR );
//    Mat m = imread("C:/Users/Laurent.PC-LAURENT-VISI/Downloads/tun.jpg",IMREAD_COLOR);
//    Mat m = imread("C:/Users/Laurent.PC-LAURENT-VISI/Downloads/5DnwY.jpg",IMREAD_COLOR);
//    Mat m = imread("C:/Users/Laurent.PC-LAURENT-VISI/Downloads/shgt_highway1.png",IMREAD_COLOR);
    
    

    double minVal,maxVal;
    float gainMax=2,gainMin=1;
    int begSlope=2*m.cols/8;
    int endSlope=3*m.cols/8;
    imshow("original", m);
    waitKey(20);
    Mat mYuv,mf;

    if (m.channels()>=3)
        cvtColor(m,mYuv,COLOR_BGR2YCrCb);
    else
        m.copyTo(mYuv);
    vector<Mat> plan;

    split(mYuv,plan);
    for (int i = 0; i < 3; i++)
    {
        minMaxLoc(plan[i],&minVal,&maxVal);
        cout << minVal << "\t"<<maxVal<<endl;
    }

    plan[0].convertTo(mf,CV_32FC1);
    mf = mf+0.01;
    log(mf,mf);

    Mat tfMf,filtre=Mat::zeros(mf.rows,mf.cols,CV_32FC1),tfMfFilter,mfHomomorphic;

    dft(mf,tfMf,cv::DFT_COMPLEX_OUTPUT);
    tfMfFilter=tfMf.clone();
    int middleWidth = mf.cols / 2 ,middleHeight = mf.rows / 2;
    double a=0.15; //Butterworth equations for homomorphic ltering of images Computers in Biology and Medicine 28 (1998) 169±181
    double d=1.5,e=0.5;
    double n=1;
    char key = 'r';
    while (key!=27)
    {
        bool newFilter=true;
        key = waitKey();
        switch (key){
        case 'd':
            d-=0.1;
            if (d==0)
                d=0.1;
            break;
        case 'D':
            d+=0.1;
            if (d==0)
                d=0.1;
            break;
        case 'e':
            e-=0.1;
            if (e==0)
                e= -0.1;
            break;
        case 'E':
            e+=0.1;
            if (e==0)
                e=0.1;
            break;
        case 'a':
            a-=0.02;
            if (a==0)
                a=-0.02;
            break;
        case 'A':
            a+=0.02;
            if (a==0)
                a=+0.02;
            break;
        case 'n':
            n -=0.1;
            if (n==0)
                n=-0.1;
            break;
        case 'N':
            n +=0.1;
            if (n==0)
                n=0.1;
            break;
        default :
            newFilter=false;
            break;
        }
        if (newFilter)
        {
            cout << "d = " << d << "\e = " << e << "\t a="<<a<< "\t n="<<n<<endl;
            for (int i = 0; i <= middleHeight; i++)
            {
                for (int j = 0; j <= middleWidth; j++)
                {
                    double r=sqrt(i*i*1.0/mf.rows*1.0/mf.rows+j*j*1.0/mf.cols*1.0/mf.cols);
                    double fLow = 1. / (1 + pow(r/a,n));
                    double fHigh =1-fLow;
                    double fBoost=d*fHigh+e;
                    if (i==0 && j==0)
                        fBoost=1;
                    filtre.at<float>(i,j)=fBoost;
                    if (i!=0)
                        filtre.at<float>(mf.rows-i,j)=filtre.at<float>(i,j);
                    if (j!=0)
                        filtre.at<float>(i,mf.cols-j)=filtre.at<float>(i,j);
                    if (j!=0 && i!=0)
                        filtre.at<float>(mf.rows-i,mf.cols-j)=filtre.at<float>(i,j);
                }
            }
    
            for (int i = 0; i <mf.rows; i++)
                for (int j = 0; j < mf.cols; j++)
                {
                    Vec2f v=tfMf.at<cv::Vec2f>(i,j);
                    complex<float> z(v[0],v[1]);
                    z = z*filtre.at<float>(i,j);
                    v[0] = z.real(); v[1]=z.imag();
                    tfMfFilter.at<Vec2f>(i,j)= v;
                }
            idft(tfMfFilter,mfHomomorphic,cv::DFT_SCALE|DFT_REAL_OUTPUT);
            minMaxLoc(mfHomomorphic,&minVal,&maxVal);
            exp(mfHomomorphic,mfHomomorphic);
            mfHomomorphic=mfHomomorphic;
            minMaxLoc(mfHomomorphic,&minVal,&maxVal);
            Mat result,r;
        //    mfHomomorphic.convertTo(r, CV_8UC1, 255 / (maxVal - minVal), 255*minVal / (maxVal-minVal));
            mfHomomorphic.convertTo(r, CV_8UC1);
            plan[0]=r;
            if (m.channels()>=3)
                merge(plan,result);
            else
                result=r;
            FileStorage fs1("homomorphic.yml", cv::FileStorage::WRITE);
	        fs1<<"Image"<<mfHomomorphic;
            cvtColor(result,r,COLOR_YCrCb2BGR);
            imshow("result",r);
            imwrite("result.jpg",r);
        }

    }


    return 0;
}
