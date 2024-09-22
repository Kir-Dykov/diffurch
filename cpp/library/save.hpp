#pragma once

#include <vector>
#include <random>
#include <string>
#include <fstream>
#include <iostream>
#include <iterator> 
#include <algorithm>

using namespace std;



template <typename T> void save(vector<T> image, string filename) {
    
    
    ofstream file(filename, ios::binary); // open file to write to
    
    // cout << "size of T is " << sizeof(T) << endl;
    if (!file) {
        cout << "Failed to open the file: " << filename << endl;
        return;
    }
    
    file.write((char*)(&(typeid(T).name()[0])), sizeof(char)); // type of array
    
    int shape_size = 1;
    int       size1= image.size();
    
	file.write((char*)(&shape_size), sizeof(int)); // dimentionality
    
	file.write((char*)(&size1), sizeof(int)); // size
    
    file.write((char*)(&image[0]), sizeof(T) * size1);    // array itself
}

template <typename T> void save(vector<vector<T>> image, string filename) {
    
    ofstream file(filename, ios::binary); // open file to write to
     if (!file) {
        cout << "Failed to open the file: " << filename << endl;
        return;
    }
    
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

template <typename T> void save(vector<vector<vector<T>>> image, string filename) {
    ofstream file(filename, ios::binary); // open file to write to
     if (!file) {
        cout << "Failed to open the file: " << filename << endl;
        return;
    }
    
    file.write((char*)(&(typeid(T).name()[0])), sizeof(char)); // type of array
    
    
    int shape_size  = 3;
    int       size1 = image.size();
    int       size2 = image[0].size();
    int       size3 = image[0][0].size();
    
	file.write((char*)(&shape_size), sizeof(int)); // dimentionality
    
	file.write((char*)(&size1), sizeof(int)); // size
	file.write((char*)(&size2), sizeof(int)); // size
	file.write((char*)(&size3), sizeof(int)); // size
    
    for (int i = 0; i < size1; i++) {
        for (int j = 0; j < size2; j++)
            file.write((char*)(&(image[i][j][0])), sizeof(T) * size3);    // array itself
    }
}
