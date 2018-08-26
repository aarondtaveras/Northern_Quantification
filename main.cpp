#include <iostream>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <streambuf>
#include <vector>
#include <fstream>
#include <iomanip>
#include <QApplication>
#include <QCharts>
#include "mainwindow.h"


using namespace std;

struct Data {
	string probe;
	string genotype;
	string strip;
	string band;
	float mean_over_area;
	vector<float> band_over_GAPDH;
	Data(string _probe, string _geno, string _strip,string _band,float _moa): probe(_probe),
	genotype(_geno),band(_band),strip(_strip),mean_over_area(_moa){};
};

struct LoadFactor {
	string genotype;
	string strip;
	float mean_over_area;
	LoadFactor(string _geno,string _strip,float _moa): genotype(_geno),strip(_strip),mean_over_area(_moa){};
};

float mean(vector<float> population){
	float sum = 0;
	for(int i=0;i<population.size();i++){
		sum += population[i];
	}
	return sum / population.size();
}

float standard_deviation(vector<float> population, float mean){
	float sum = 0;
	for(int i=0;i<population.size();i++){
		sum += pow((population[i] - mean), 2.0);
	}
	return sqrt(sum/population.size());
}

float tscore(vector<float> population, float s_mean, float p_mean){
	
	float s_dev = standard_deviation(population,s_mean);

	return (s_mean - p_mean) / (s_dev / sqrt(population.size()/2) );
}

void calc_band_over_GAPDH(vector<Data> & trials, vector<LoadFactor> & load_factors){

	vector<float> result;
		// Band/GAPDH
		for(auto & data: trials)
		{
			for(auto & load: load_factors)
			{
				if(data.genotype == load.genotype && data.strip == load.strip)
				{
					// band_over_GAPDH += to_string(data.mean_over_area / load.mean_over_area) + ",";
					data.band_over_GAPDH.push_back(data.mean_over_area / load.mean_over_area);
				}
			}
		}

}


void export_to_csv(ostream & o,vector<Data> trials){
	// Gotta export , vector<float> tscores, vector<float> t_mean as well! 
	for(auto & item: trials)

	{
		o << item.probe << " " << item.genotype << " ";
		for(auto & num: item.band_over_GAPDH){
		o << to_string(num) << "," << endl;	
		} 
	}
}

int main(int argc, char * argv [])
{

	ifstream input;
	ofstream output;

	input.open(argv[1]);
	if(input.fail()){
		cerr << "Could not open supplied input file. Try again";
	}

	output.open("results.txt");
	if(output.fail()){
		cerr << "Could not open specified output stream.";
	}

	vector<LoadFactor> load_factors;
	vector<Data> trials;

	string line;

	/* Parsing the CSV file. */
	while(getline(input, line)){
		stringstream trial_stream(line);
		vector<string> field;
		while( trial_stream.good()){
		string sub;
			while(getline(trial_stream,sub,',')){
				sub.erase(remove(sub.begin(),sub.end(),' '),sub.end() );
				field.push_back(sub);
				// cout << sub << "\t";
			}
			if(field[0]=="ITS1A" || field[0]=="ITS2B")
			{
				trials.push_back(Data(field[0],field[1],field[2],field[3],atof(field[9].c_str() ) ) );
			}
			else if(field[0]=="GAPDH")
			{
				load_factors.push_back(LoadFactor( field[1],field[2],atof(field[9].c_str() ) ) );
			}
			// cout << endl;
		}
	}

	/* Beginning our calculations. 
		- First is band/GAPDH
		- Second is t-score test per band. Population = all bands per type
		- Third is average of all the T-scores for given band
	*/
	if(trials.size()!=0 && load_factors.size()!=0)
	{
		// Band/GAPDH done!
		calc_band_over_GAPDH(trials,load_factors);

		// exporting to CSV ...
		export_to_csv(output, trials);
		//T-score test, but the population needs to be all of same band
		//and the sample needs to be of the mutants. 

		// for(auto & data: trials){
		// 	if(data)
		// }


	}

	
	// for(auto & d: trials){
	// 	cout << d.probe << endl; 
	// }
	QApplication a(argc, argv);
    MainWindow w;
    w.show();
    
    return a.exec();


	input.close();
	output.close();

	return 0;
}