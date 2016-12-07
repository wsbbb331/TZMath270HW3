//
//  ThreeDOutput.h
//  Math270AHW2
//
//  Created by Tianyang Zhang on 12/4/16.
//  Copyright © 2016 UCLA. All rights reserved.
//

#ifndef ThreeDOutput_h
#define ThreeDOutput_h

#include <Eigen/Dense>
#include <Eigen/Core>
#include <iostream>


namespace JIXIE {
    template <class T>
    
    class Cube {
        //unit cube class
    public:
        T scale = 1.0;
        std::vector<std::vector<T>> vertices;
        std::vector<std::vector<int>> faces;
        T l, w, h;
        inline void initCubeValue(){
            vertices = {{0.0, 0.0, 0.0},
                {0.0, 0.0, h * scale},
                {0.0, w * scale, 0.0},
                {0.0, w * scale, h * scale},
                {l * scale, 0.0, 0.0},
                {l * scale, 0.0, h * scale},
                {l * scale, w * scale, 0.0},
                {l * scale, w * scale, h * scale}};
            faces = {{1, 7, 5},
                {1, 3, 7},
                {1, 4, 3},
                {1, 2, 4},
                {3, 8, 7},
                {3, 4, 8},
                {5, 7, 8},
                {5, 8, 6},
                {1, 5, 6},
                {1, 6, 2},
                {2, 6, 8},
                {2, 8, 4}};
        }
        inline Cube(const T x,const T y,const T z, const T scaleInput = 1.0, const T lInput = 1.0, const T winput = 1.0, const T hinput = 1.0): scale(scaleInput), l(lInput), w(winput), h(hinput){
            initCubeValue();
            for (int i = 0; i < 8; i++){
                vertices[i][0] += (x - 0.5 * scale * l);
                vertices[i][1] += (y - 0.5 * scale * w);
                vertices[i][2] += (z - 0.5 * scale * h);
            }
        }
    };
    
    template <class T>
    class ObjBody {
        typedef Eigen::Matrix<T,Eigen::Dynamic,1> TVect;
    public:
        std::vector<std::vector<T>> verticesList;
        std::vector<std::vector<int>> faceList;
        inline ObjBody(){}
        inline void addCube(Cube<T>& cube){
            std::vector<T> cubeVericesIndex;
            std::vector<std::vector<int>> cubeFaceOutput = cube.faces;
            for (int i = 0; i < 8; i++){
                bool duplicateVertex = false;
                std::vector<T> newVertices = cube.vertices[i];
                std::vector<T> newVertices2 = newVertices;
                std::vector<T> newVertices3 = newVertices;
                std::rotate(newVertices2.begin(),newVertices2.begin()+1,newVertices2.end());
                for (int j = 0; j < verticesList.size(); j++){
                    if (newVertices == verticesList[j]){
                        mapNewVertexIndex(i+1, j+1, cube.faces, cubeFaceOutput);
                        duplicateVertex = true;
                        break;
                    }
                }
                if(!duplicateVertex){
                    verticesList.push_back(newVertices);
                    mapNewVertexIndex(i+1, int(verticesList.size()), cube.faces, cubeFaceOutput);
                }
            }
            
            for(int i=0;i<12;i++){
                bool duplicateFace = false;
                std::vector<int> newFace = cubeFaceOutput[i];
                for (int j=0; j < faceList.size(); j++){
                    if (is_permutation(newFace.begin(), newFace.end(), faceList[j].begin())){
                        duplicateFace = true;
                        break;
                    }
                }
                if(!duplicateFace){
                    faceList.push_back(newFace);
                }
            }
        }
        
        inline void mapNewVertexIndex(const int& oldIndex, const int& newIndex, const std::vector<std::vector<int>> faceSearchList, std::vector<std::vector<int>>& faceList){
            for (int i=0; i<faceSearchList.size(); i++){
                for (int j=0; j<faceSearchList[i].size(); j++){
                    if(faceSearchList[i][j] == oldIndex)(faceList[i][j] = newIndex);
                }
            }
        }
        
        inline void printOutput(){
            for (int i=0;i<verticesList.size();i++){
                std::cout<<"v ";
                for (int j=0;j<verticesList[i].size();j++){
                    std::cout << " "<< verticesList[i][j];
                }
                std::cout << std::endl;
            }
            for (int i=0;i<faceList.size();i++){
                std::cout<<"f ";
                for (int j=0;j<faceList[i].size();j++){
                    std::cout<< " "<<faceList[i][j];
                }
                std::cout<<std::endl;
            }
        }
        
        inline void writeOutput(std::string filename){
            std::ofstream outdata;
            std::string simulation_data_filename(std::string("output/obj/")+filename+std::string(".obj"));
            outdata.open(simulation_data_filename.c_str());
            for (int i=0;i<verticesList.size();i++){
                outdata<<"v ";
                for (int j=0;j<verticesList[i].size();j++){
                    outdata << " "<< verticesList[i][j];
                }
                outdata << std::endl;
            }
            for (int i=0;i<faceList.size();i++){
                outdata<<"f ";
                for (int j=0;j<faceList[i].size();j++){
                    outdata<< " "<<faceList[i][j];
                }
                outdata<<std::endl;
            }
            outdata.close();
        }
    
    };
    
    template <class T>
    class Shape {
    public:
        std::vector<Cube<T>> cubeVector;
        inline Shape(const T l, const T w, const T h, const T scale, const T x=0.0, const T y=0.0, const T z=0.0){
            for (T i=-l/2.0+x+0.5*scale;i<=l/2.0+x-0.5*scale;i=i+scale){
                for (T j=-w/2.0+y+0.5*scale;j<=w/2.0+y-0.5*scale;j=j+scale){
                    for (T k=-h/2.0+z+0.5*scale;k<=h/2.0+z-0.5*scale;k=k+scale){
//                        std::cout<<i<<" "<<j<<" "<<k<<std::endl;
                        Cube<T> cubeElement(i,j,k,scale);
                        cubeVector.push_back(cubeElement);
                    }
                }
            }
        }
        inline void insertObj(ObjBody<T>& body){
            for (int i=0; i<cubeVector.size(); i++){
                body.addCube(cubeVector[i]);
            }
        }
    };
    
    
}



#endif /* ThreeDOutput_h */
