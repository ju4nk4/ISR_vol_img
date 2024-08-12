// 2022-05-18, Juan Carlos Araújo, ju4nk4@gmail.com

/*

H5::Group
H5::DataSpace
H5::DataSet
H5::PredType

*/

// https://stackoverflow.com/questions/3746484/why-am-i-getting-this-redefinition-of-class-error


#include <vector>
#include <iostream>
#include <memory>
#include "H5Cpp.h"
//#include "H5Exception.h"


// #pragma once

struct Dim {
    hsize_t cols, rows;
};

Dim getDimensions(H5::DataSet& set);
	


struct BeamInfo {
	unsigned int id, index_i, index_f, num_range;
	double elm, azm;
};


template<class Type>
class Elem2D {
private:
	Type* _data;
	int _base;

public:
    Elem2D(int base, Type* data) : _data(data), _base(base) {}
	
    Type operator[] (int index) {
        return _data[_base + index];
    }
    void clear () {
    	delete[] _data;
    }
};


template<class Type>
class Array2D {
public:
    Array2D(){}
    Array2D(Type* data, int size) : _data(data),_size(size) {}

    Elem2D<Type> operator[] (int index) {
        return Elem2D<Type>(index*_size, _data);
    }
    void clear () {
    	delete[] _data;
    }
private:
    Type* _data;
    int _size;
};


//----------------------------------------------------------------------------------------
class Metadata {
	public:
		std::string instrument;
		std::string inst_code;
		std::string inst_kind;
		std::string start_t;
		std::string end_t;
		std::string cedar_name;
		std::string inst_lat;
		std::string inst_lon;
		std::string inst_alt;
		std::string inst_cat;
		
		void print () {
		
			printf ("Instrument: \t%s\n", instrument.c_str() );
			printf ("Code: \t%s\n", inst_code.c_str() );
			printf ("Kind: \t%s\n", inst_kind.c_str() );
			printf ("Start time: \t%s\n", start_t.c_str() );
			printf ("End time: \t%s\n", end_t.c_str() );
			printf ("Cedar: \t%s\n", cedar_name.c_str() );
			printf ("Inst latitude: \t%s\n", inst_lat.c_str() );
			printf ("Inst longitude: \t%s\n", inst_lon.c_str() );
			printf ("Inst altitude: \t%s\n", inst_alt.c_str() );
			printf ("Category: \t%s\n", inst_cat.c_str() );
		}
};

void read_metadata ( std::string path, Metadata &meta );


//template<class NativeType>
		//NativeType* 
void read1DStrData(const H5::PredType& type, const H5::DataSet data, hsize_t offset, hsize_t offset_size, hsize_t size);


template<class NativeType>
NativeType* read1DData(const H5::PredType& type, const H5::DataSet data, hsize_t offset, hsize_t offset_size, hsize_t size) {

    hsize_t     dimsm[1];
    dimsm[0] = size;
    H5::DataSpace memspace(1, dimsm);

    hsize_t      offset_out[1];   // hyperslab offset in memory
    hsize_t      count_out[1];    // size of the hyperslab in memory
    offset_out[0] = offset;
    count_out[0] = offset_size;
    memspace.selectHyperslab(H5S_SELECT_SET, count_out, offset_out);

    NativeType* out_data = new NativeType[size];
    H5::DataSpace space = data.getSpace();
    data.read(out_data, type, memspace, space);
    return out_data;

}


template<class NativeType>
Array2D<NativeType> read2DData(const H5::PredType& type, const H5::DataSet& data, const Dim& offset, const Dim& offset_size, const Dim& size) {
    H5::DataSpace space = data.getSpace();
    hsize_t      offset_arr[2];   // hyperslab offset in the file
    hsize_t      count_arr[2];    // size of the hyperslab in the file
    offset_arr[0] = offset.cols;
    offset_arr[1] = offset.rows;
    count_arr[0] = offset_size.cols;
    count_arr[1] = offset_size.rows;
    space.selectHyperslab(H5S_SELECT_SET, count_arr, offset_arr);

    hsize_t     dimsm[2];
    dimsm[0] = size.cols;
    dimsm[1] = size.rows;
    H5::DataSpace memspace(2, dimsm);

    hsize_t      offset_out[2];   // hyperslab offset in memory
    hsize_t      count_out[2];    // size of the hyperslab in memory
    offset_out[0] = offset.cols;
    offset_out[1] = offset.rows;
    count_out[0] = offset_size.cols;
    count_out[1] = offset_size.rows;
    memspace.selectHyperslab(H5S_SELECT_SET, count_out, offset_out);

    NativeType *out_data = new NativeType[size.cols*size.rows];
    try {
        data.read(out_data, type, memspace, space);
    } catch(const std::exception& err) {
        exit(0);
    }
    return Array2D<NativeType>(out_data,size.rows);

}

std::vector<H5::Group> getAllGroups(H5::Group& group);




