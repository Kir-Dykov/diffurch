#pragma once

#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator> 
#include <algorithm>

using namespace std;

template <typename T> void save(vector<vector<T>> image, string filename) {
    ofstream file(filename, ios::binary); // open file to write to
    file.write((char*)(&(typeid(T).name()[0])), sizeof(char)); // type of array
    
    
    int shape_size  = 2;
    int       size1 = image.size();
    int       size2 = image[0].size();
    
	file.write((char*)(&shape_size), sizeof(int)); // dimentionality
    
	file.write((char*)(&size1), sizeof(int)); // size
	file.write((char*)(&size2), sizeof(int)); // size
    
    for (int i = 0; i < size1; i++) {
        file.write((char*)(&(image[i][0])), sizeof(T) * size2);    // array itself
    }
}

template <typename T> void save(vector<T> image, string filename) {
    ofstream file(filename, ios::binary); // open file to write to
    file.write((char*)(&(typeid(T).name()[0])), sizeof(char)); // type of array
    
    // cout << "size of T is " << sizeof(T) << endl;
    
    int shape_size = 1;
    int       size1= image.size();
    
	file.write((char*)(&shape_size), sizeof(int)); // dimentionality
    
	file.write((char*)(&size1), sizeof(int)); // size
    
    file.write((char*)(&image[0]), sizeof(T) * size1);    // array itself
}


