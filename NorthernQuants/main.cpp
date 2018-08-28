#include <iostream>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <streambuf>
#include <vector>
#include <fstream>
#include <iomanip>
#include <QApplication>
#include <QtCharts>
#include "mainwindow.h"

using namespace std;

struct Strip {
    string probe;
    bool is_load_control;
    string genotype;
    string id;
    string band;
    float mean_over_area;
    vector<float> band_over_GAPDH;
    Strip(string _probe, string _geno, string _id,string _band,float _moa): probe(_probe),
    genotype(_geno),band(_band),id(_id),mean_over_area(_moa){}
};

struct LoadingControl {
    string genotype;
    string id;
    float mean_over_area;
    LoadingControl(string _geno,string _id,float _moa): genotype(_geno),id(_id),mean_over_area(_moa){}
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

void calc_band_over_GAPDH(vector<Strip> & strips, vector<LoadingControl> & loading_controls){

    vector<float> result;
        // Band/GAPDH
        for(auto & data: strips)
        {
            for(auto & load: loading_controls)
            {
                if(data.genotype == load.genotype && data.id == load.id)
                {
                    // band_over_GAPDH += to_string(data.mean_over_area / load.mean_over_area) + ",";
                    data.band_over_GAPDH.push_back(data.mean_over_area / load.mean_over_area);
                }
            }
        }

}


void export_to_csv(ostream & o,vector<Strip> strips){
    // Gotta export , vector<float> tscores, vector<float> t_mean as well! 
    for(auto & item: strips)

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

    vector<LoadingControl> loading_controls;
    vector<Strip> strips;

    string line;

    /* Parsing the CSV file. */
    while(getline(input, line)){
        stringstream strip_stream(line);
        vector<string> field;
        while( strip_stream.good()){
        string sub;
            while(getline(strip_stream,sub,',')){
                sub.erase(remove(sub.begin(),sub.end(),' '),sub.end() );
                field.push_back(sub);
                // cout << sub << "\t";
            }
            if(field[0]=="ITS1A" || field[0]=="ITS2B")
            {
                strips.push_back(Strip(field[0],field[1],field[2],field[3],atof(field[9].c_str() ) ) );
            }
            else if(field[0]=="GAPDH")
            {
                loading_controls.push_back(LoadingControl( field[1],field[2],atof(field[9].c_str() ) ) );
            }
            // cout << endl;
        }
    }

    /* Beginning our calculations. 
        - First is band/GAPDH
        - Second is t-score test per band. Population = all bands per type
        - Third is average of all the T-scores for given band
    */
    if(strips.size()!=0 && loading_controls.size()!=0)
    {
        // Band/GAPDH done!
        calc_band_over_GAPDH(strips,loading_controls);

        // exporting to CSV ...
        export_to_csv(output, strips);
        //T-score test, but the population needs to be all of same band
        //and the sample needs to be of the mutants. 

        // for(auto & data: strips){
        //  if(data)
        // }


    }

    input.close();
    output.close();

    QApplication a(argc, argv);
    QBarSet *set0 = new QBarSet("ITS1A");
    QBarSet *set1= new QBarSet("ITS2B");
    QBarSet *LOAD= new QBarSet("LOAD");
    QBarSet *GAPDH = new QBarSet("GAPDH");

    QStackedBarSeries *series = new QStackedBarSeries();
    series->append(set0);
    series->append(set1);

    QChart *chart = new QChart();
    chart -> addSeries(series);
    chart -> setTitle("Tester");
    chart -> setAnimationOptions(QChart::SeriesAnimations);

    QStringList categories;
    categories << "Jan" << "Feb" << "Mar" << "Apr" << "May" << "Jun";
    QBarCategoryAxis *axis = new QBarCategoryAxis();
    axis->append(categories);
    chart->createDefaultAxes();
    chart->setAxisX(axis, series);

    chart->legend()->setVisible(true);
    chart->legend()->setAlignment(Qt::AlignBottom);

    QChartView *chartView = new QChartView(chart);
    chartView->setRenderHint(QPainter::Antialiasing);

    QMainWindow w;
    w.setCentralWidget(chartView);
    w.show();

    return a.exec();

}
