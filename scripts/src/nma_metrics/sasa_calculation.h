# ifndef _SASA_CALCULATION_H
# define _SASA_CALCULATION_H

# include<map>
# include<cmath>
# include<vector>
# include<cassert>
# include<cstring>
# include<sstream>
# include<iostream>
# include<dlib/dnn.h>
# include<dlib/data_io.h>

# include "amino.h"
# include "utility.h"
# include "geometry.h"
# include "trajectory.h"

using namespace std;
using namespace dlib;

# ifndef INPUT_LINE_LENGTH
# define INPUT_LINE_LENGTH  1000
# endif

# ifndef MIN_SAMPLES
# define MIN_SAMPLES 10
# endif

# ifndef NNBRS
# define NNBRS 15
# endif

struct RegressorInput{
public:
	std::vector<matrix<float>> x;
	std::vector<float> y;
public:
	RegressorInput():x(),y()
	{}
	
	RegressorInput& operator = (RegressorInput const& rhs)
	{
		this->x  = rhs.x;
		this->y  = rhs.y;
		return *this;
	}
	
	int ndata()const
	{ 
		return static_cast<int>(this->x.size());
	}
};

using sasanet = loss_mean_squared<fc<1,relu<fc<15,relu<fc<60,relu<fc<60,input<matrix<float>>>>>>>>>>;

RegressorInput read_dnn_regressor_input( string const& );

void train_test_split(RegressorInput const&, RegressorInput&, RegressorInput&, float = 0.2);

float train_network(RegressorInput const& , string const& , float = 0.2);


RegressorInput read_dnn_regressor_input( string const& filename )
{
	assert( ::is_file(filename.c_str()) );
	
	char   line[ INPUT_LINE_LENGTH + 1];
	bzero( line, INPUT_LINE_LENGTH + 1);
	int m, n = -1;
	ifstream f(filename);
    assert( f.good() );
	std::vector< std::vector<float> > data;
	while( f.getline(line, INPUT_LINE_LENGTH) ){
		std::vector<string> words = string_split(string(line),",");
		if( n < 0) n = words.size();
		assert( n > 1 && n == words.size() );
		std::vector<float>  row;
		for(int i=0; i < words.size(); ++i){
			::is_float( words[i] );
			row.push_back( atof(words[i].c_str()) );
		}
		data.push_back(row);
		bzero( line, INPUT_LINE_LENGTH );
	}
	f.close();
	m = data.size();
	RegressorInput inp;
	for(int i=0; i < m; ++i )
	{
		inp.x.push_back(dlib::matrix<float>(n-1,1));
		inp.y.push_back(data[i][n-1]);
		for(int j=0; j < n-1; ++j )
			inp.x[i](j,0) = data[i][j];
		
	}
	return inp;
}

void train_test_split(RegressorInput const& data, 
                      RegressorInput& train, 
					  RegressorInput& test, 
					  float testSplit)
{
	assert( testSplit < 0.9 && testSplit >= 0);
	int n = data.ndata();
	assert( n > MIN_SAMPLES );
	int test_size = static_cast<int>( n * testSplit );
	int train_size = n - test_size;
	std::vector<int> test_index, train_index;
	int pass = 0;
	for( int i=0; test_index.size() < test_size; i++ ){
		if( train_index.size() < train_size && test_index.size() < test_size){
			if( ::rand()/(RAND_MAX + 1.f) < testSplit )
				test_index.push_back(i);
			else
				train_index.push_back(i);
		}else if( train_index.size() < train_size ){
			train_index.push_back(i);
		}else{
			test_index.push_back(i);
		}
	}
	
	train.x.clear(); train.y.clear();
	test.x.clear();  test.y.clear();
	
	for(int i=0; i < test_index.size(); ++i){
		test.x.push_back( data.x[test_index[i]] );
		test.y.push_back( data.y[test_index[i]] );
	}
	
	for(int i=0; i < train_index.size(); ++i){
		train.x.push_back(data.x[train_index[i]]);
		train.y.push_back(data.y[train_index[i]]);
	}
}


float train_network( RegressorInput const& inp, string const& model_file, float test_split )
{
	sasanet nn;
	adam opt;
	
	RegressorInput test_input, train_input;
	
	train_test_split(inp, train_input, test_input, test_split);
	
	dnn_trainer<sasanet,adam> trainer(nn,opt);
	trainer.set_learning_rate(0.01);
	trainer.set_min_learning_rate(1e-5);
	trainer.set_mini_batch_size(1000);
	trainer.be_verbose();
	trainer.train(train_input.x, train_input.y);
	nn.clean();
	serialize(model_file.c_str()) << nn;
	
	float err = 0.f;
	int n = test_input.ndata();
	
	std::vector<float> predicted_value = nn(test_input.x);
	assert( predicted_value.size() == test_input.ndata() );
	for(int i=0; i < n; ++i )
		err += fabs( predicted_value[i] - test_input.y[i]);
	return err/n ;
}


class SASAPredictor
{
protected:
	static SASAPredictor*  m_inst;
	std::map<string, sasanet> m_models;

private:
    SASAPredictor():m_models()
	{
		string model_location = ::base_directory() + "/data/";
		int n = sizeof(Amino3Letter)/sizeof(string);
		for( int i = 0; i < n; ++i )
		{
			std::ostringstream  fname;
			fname.clear();
			fname << Amino3Letter[i] <<".dat";
			string model_file = model_location + fname.str();
			if( ::is_file(model_file) )
			{
			   sasanet nn;
			   deserialize(model_file.c_str()) >> nn;
			   this->m_models[Amino3Letter[i]] = nn;
			}
		}
	}

public:
    static SASAPredictor& get_instance()
	{
		if( SASAPredictor::m_inst == nullptr )
			SASAPredictor::m_inst = new SASAPredictor();
		return *m_inst;
	}
	
	float predict( string const& aa, dlib::matrix<float> const& feature )
	{
		assert( ::is_amino( aa ) );
		assert( ::is_in(this->m_models, aa) );
		std::vector<dlib::matrix<float> > inp;
		inp.push_back( feature );
		return this->m_models[aa](inp)[0];
	}
	
	std::vector<float> predict( string const& aa, std::vector<dlib::matrix<float>> const& features )
	{
		assert( ::is_amino( aa ) );
		assert( ::is_in(this->m_models, aa ) );
		return this->m_models[aa](features);
	}

};

SASAPredictor* SASAPredictor::m_inst = nullptr;


std::map<int, dlib::matrix<float> > snapshot_to_sasa_feature( Snapshot const& snap, int topn = NNBRS )
{
	int n = ::snapshot_size( snap );
	assert( n > 2 );
	std::vector<int> resids = snap.residueIds;
	std::map<int, std::vector<NeighborData> > nbrs = ::top_n_neighbors(snap, topn);
	std::map<int, dlib::matrix<float> > res;
    for( int i=1; i < n-1; ++i )
	{
		Vector3D du = (::connecting_vector(snap.coordinates[i-1], snap.coordinates[i]).unit_vector() + 
		              ::connecting_vector(snap.coordinates[i+1], snap.coordinates[i]).unit_vector()).unit_vector();
		assert( ::is_in( nbrs, resids[i] ) );
		std::vector<NeighborData> ndata = nbrs[resids[i]];
		assert( ndata.size() == topn );
		std::vector<float> vfeatures;
		for( int j=0; j < topn; ++j)
		{
			Vector3D dv = ::connecting_vector(snap.coordinates[i], snap.coordinates[ndata[j].nbr_residx] ).unit_vector();
			float c = ::dot(dv,du);
			float d = ndata[j].nbr_distance;
			float v = ::get_volume(snap.residueType[ndata[j].nbr_residx]);
			float h = ::get_hydropathy(snap.residueType[ndata[j].nbr_residx]);
			vfeatures.push_back(c);
			vfeatures.push_back(h); //Change this line after fixing the bug
			vfeatures.push_back(h);
			vfeatures.push_back(d);
		}
		res[ resids[i] ] = dlib::matrix<float>(vfeatures.size(),1);
		for(int j=0; j < vfeatures.size(); ++j )
			res[ resids[i] ](j,0) = vfeatures[j];
	}	
	return res;
}

std::map<int,float> sasa_values( Snapshot const& snapshot, std::vector<int> const& resids )
{
	assert( resids.size() > 0 );
	std::map<int, dlib::matrix<float> > res = ::snapshot_to_sasa_feature(snapshot);
	SASAPredictor& predictor = SASAPredictor::get_instance();
	std::map<int,float> sasa_map;
	
	for(typename std::vector<int>::const_iterator i = resids.begin(); i != resids.end(); ++i )
	{
		string restype = ::fetch_resname(snapshot,*i);
		if( ::is_in( res, *i) )
		{
			sasa_map[*i] = predictor.predict(restype, res[*i]);
		}
	}
	return sasa_map;
}
	

# endif