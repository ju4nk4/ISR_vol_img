// 2022-05-18, Juan Carlos Araújo, ju4nk4@gmail.com

#include "H5Cpp.h"
//#include "H5Exception.h"

using namespace H5;

struct Dim {
    hsize_t cols, rows;
};

Dim getDimensions(DataSet& set) {
    DataSpace space = set.getSpace();
    int rank = space.getSimpleExtentNdims();
    if (rank == 1) {
        hsize_t dims[1];
        space.getSimpleExtentDims(dims, NULL);
        return { dims[0],1 };
    }
    hsize_t dims[2];
    space.getSimpleExtentDims(dims,NULL);
    return { dims[0],dims[1] };
}
	


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

void read_metadata ( string path, Metadata &meta ) {
// https://stackoverflow.com/questions/20778978/difficulty-parsing-hdf5-compound-data-type
// https://stackoverflow.com/questions/13814027/reading-a-string-from-hdf5-in-c


	/* // This is the structure of the H5T_COMPOUND
	  DATATYPE  H5T_COMPOUND {
		  H5T_STRING {
			 STRSIZE 20;
			 STRPAD H5T_STR_NULLPAD;
			 CSET H5T_CSET_ASCII;
			 CTYPE H5T_C_S1;
		  } "name";
		  H5T_STRING {
			 STRSIZE 67;
			 STRPAD H5T_STR_NULLPAD;
			 CSET H5T_CSET_ASCII;
			 CTYPE H5T_C_S1;
		  } "value";
	   }
	   DATASPACE  SIMPLE { ( 14 ) / ( 14 ) }
	*/
	
	class RawData {
		public:
			char name [100]; // STRSIZE 20;
			char value[100]; // STRSIZE 67;
	};
	
	printf ("Reading file\n");
	string ROOT_GROUP = "/Metadata", PAR_1D = "Experiment Parameters";
	
	/*
	std::string 
		directory = "/home/juan/postdoc-visualization-files/ionosradar_data/data",
		filename  = "ras161121.002",
		path = directory + "/" + filename + ".hdf5";
	*/
	
	printf ("path: %s\n", path.c_str() );
	std::unique_ptr<H5File>     _file = nullptr;
	const H5std_string h5filename(path);
    try {
        _file = std::make_unique<H5File>(path.c_str(), H5F_ACC_RDONLY);
    } catch (const H5::FileIException& err) {
		std::cerr << "Info: Error reading HDF5 data: {} " << std::endl;
    }
	
	Group data                  = _file->openGroup(ROOT_GROUP);
	
	auto dataset = data.openDataSet( PAR_1D );
	hsize_t numOfmetadata       = getDimensions(dataset).cols; 
	//printf ("***** %lld ******\n\n", numOfmetadata );
    
    auto dataspace = dataset.getSpace();
    //hsize_t dims_out[2];
    // auto ndims = dataspace.getSimpleExtentDims(dims_out, nullptr);
    
    // auto n = dims_out[0] * dims_out[1];
    auto data_type = dataset.getDataType();
    //auto type_class = data_type.getClass();
    //auto data_size = data_type.getSize();
	
	
	
	CompType type(sizeof(RawData));
	type.insertMember("name",  HOFFSET(RawData, name ), StrType(0, 100));
	type.insertMember("value", HOFFSET(RawData, value), StrType(0, 100));
	
	
	std::vector<RawData> data_vector;
	data_vector.resize( numOfmetadata );  // 14
	
	dataset.read( data_vector.data(), type);
	
	unsigned int counter = 0;
	for( RawData &m : data_vector) {
		// printf ("counter = %d, (%s, %s) \n", counter, m.name, m.value );
		//printf ("\t(%s, %s) \n", m.name, m.value );
		
		
		if (strcmp(m.name, "instrument" ) == 0) {
		    meta.instrument = string(m.value);
		    //std::cout << "\t\tFound instrument = " << meta.instrument << std::endl;
		}
		if (strcmp(m.name, "instrument code(s)" ) == 0) {
		    //std::cout << "\t\tFound instrument code(s) = " << m.value << std::endl;
		    meta.inst_code = string(m.value);
		}
		if (strcmp(m.name, "kind of data file" ) == 0) {
		    //std::cout << "\t\tFound kind of data file = " << m.value << std::endl;
		    meta.inst_kind = string(m.value);
		}
		if (strcmp(m.name, "instrument category" ) == 0) {
		    //std::cout << "\t\tFound instrument category = " << m.value << std::endl;
		    meta.inst_cat = string(m.value);
		}
		if (strcmp(m.name, "Cedar file name" ) == 0) {
		    //std::cout << "\t\tFound Cedar file name = " << m.value << std::endl;  
		    meta.cedar_name = string(m.value);
		}
		if (strcmp(m.name, "start time" ) == 0) {
		    //std::cout << "\t\tFound start time = " << m.value << std::endl;
		    meta.start_t = string(m.value);
		    //read_UT_time (m.value);
		}
		if (strcmp(m.name, "end time" ) == 0) {
		    //std::cout << "\t\tFound end time = " << m.value << std::endl;
		    meta.end_t = string(m.value);
		}
		if (strcmp(m.name, "instrument latitude" ) == 0) {
		    //std::cout << "\t\tFound latitude = " << m.value << std::endl;
		    meta.inst_lat = string(m.value);
		}
		if (strcmp(m.name, "instrument longitude" ) == 0) {
		    //std::cout << "\t\tFound longitude = " << m.value << std::endl;
		    meta.inst_lon = string(m.value);
		}
		if (strcmp(m.name, "instrument altitude" ) == 0) {
		    //std::cout << "\t\tFound altitude = " << m.value << std::endl;
		    meta.inst_alt = string(m.value);
		}
		
		counter++;
	}
}


//template<class NativeType>
		//NativeType* 
void read1DStrData(const PredType& type, const DataSet data, hsize_t offset, hsize_t offset_size, hsize_t size) {
	
	// ************** NOT WORKING ************************************************
	
	hsize_t     dimsm[1];
	dimsm[0] = size;
	DataSpace memspace(1, dimsm);

	hsize_t      offset_out[1];   // hyperslab offset in memory
	hsize_t      count_out[1];    // size of the hyperslab in memory
	offset_out[0] = offset;
	count_out[0] = offset_size;
	memspace.selectHyperslab(H5S_SELECT_SET, count_out, offset_out);
	
	std::vector<string> out_data(size);
	
	//NativeType* out_data = new NativeType[size];
	DataSpace space = data.getSpace();
	
	string out;
	data.read(out, type, memspace, space);
	// return out_data;
	
	// ************** NOT WORKING ************************************************
}


template<class NativeType>
NativeType* read1DData(const PredType& type, const DataSet data, hsize_t offset, hsize_t offset_size, hsize_t size) {

    hsize_t     dimsm[1];
    dimsm[0] = size;
    DataSpace memspace(1, dimsm);

    hsize_t      offset_out[1];   // hyperslab offset in memory
    hsize_t      count_out[1];    // size of the hyperslab in memory
    offset_out[0] = offset;
    count_out[0] = offset_size;
    memspace.selectHyperslab(H5S_SELECT_SET, count_out, offset_out);

    NativeType* out_data = new NativeType[size];
    DataSpace space = data.getSpace();
    data.read(out_data, type, memspace, space);
    return out_data;

}



template<class NativeType>
Array2D<NativeType> read2DData(const PredType& type, const DataSet& data, const Dim& offset, const Dim& offset_size, const Dim& size) {
    DataSpace space = data.getSpace();
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
    DataSpace memspace(2, dimsm);

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

std::vector<Group> getAllGroups(Group& group) {
    std::vector<Group> vals;
    const int chr_size = 128;
    int size = group.getNumObjs();
    for (int i = 0; i < size; i++) {
        char name[chr_size];
        int name_size = group.getObjnameByIdx(i, &name[0], chr_size);
        vals.push_back( group.openGroup(std::string(name, name_size).c_str()) );
        // printf("%%\t%s\n", std::string(name, name_size).c_str() );
    }
    return vals;
}






