#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<cmath>
#include<complex>
#include <vector>
#include <ctime>

#define PREC double

struct Site;
template<class floatT> class GridAccessor;
template<class floatT> class Grid;

void indexerTest(int max);

// structure to define the lattice points leaving out the z component due to the contemplated event plane
struct Site{
       std::vector<int> VCoordinates;

       // constructor
       Site(int x, int y){
              VCoordinates.push_back(x);
              VCoordinates.push_back(y);
       }

};

// class to store data
template<class floatT>
class Grid{
       private:
              // vector to store the data
              std::vector<floatT> _VData;

              // number of sites in the 3 spartial directions (only cubic grids are allowed)
              int _maxSitesPerDirection;
       public:
              //constructor
              Grid(int maxSitesPerDirection) :
                                   _VData( 2*maxSitesPerDirection
                                          + maxSitesPerDirection*maxSitesPerDirection
                                          + 1 ) ,
                                   _maxSitesPerDirection(maxSitesPerDirection) {}

              int getIndex(Site site){
                     int max = _maxSitesPerDirection;
                     return site.VCoordinates[0] + site.VCoordinates[1] + site.VCoordinates[1]*max;
              }

              void setSite(Site site, floatT value){
                     _VData.at(getIndex(site)) = value;
              }

              floatT getSite(Site site){
                     return _VData.at(getIndex(site));
              }

              void addSite(Site site, floatT value){
                     _VData.at(getIndex(site)) += value;
              }

              int getMaxSitesPerDirection(){
                     return _maxSitesPerDirection;
              }
};


template<class floatT>
class FileWriter{
       private:
              int _NEvents;
              int _NNucleonsCore;
              std::string _filename;
       public:
              // constructor
              FileWriter(int newNEvents, int newNNucleonsCore, std::string newFilename) : _NEvents(newNEvents),
                                                        _NNucleonsCore(newNNucleonsCore),
                                                        _filename(newFilename) {}

              void readFile(Grid<floatT> & grid){
                     std::fstream file;
                     file.open(_filename.c_str(), std::ios::in);

                     if(!file.is_open()){return;}

                     std::vector<floatT> row(4);

                     // loop through the data file which format is clarified by
                     // x corrd \t y coord \t z coord \t NColl
                     for(int i = 0; i < 2 * _NEvents * _NNucleonsCore; ++i){
                            for(int col = 0; col < 4 ; ++col){
                                   file >> row[col];
                            }
                            // translate the euclidean coordinates into the grid coordinates
                            // leaving out the z component
                            int x = (row[0] + 15)*100;
                            int y = (row[1] + 15)*100;

                            if(x > grid.getMaxSitesPerDirection() ||
                               y > grid.getMaxSitesPerDirection()  ) {
                                          std::cout << "ERROR@readFile: Coordinates out of grid!" << '\n';
                            }

                            // compute the site
                            Site site(x, y);

                            // set NColl on the grid
                            grid.setSite(site, row[3]);
                     }

                     file.close();
              }

              void writeFileEDens(Grid<floatT> * grid, int accuracy = 10){
                     std::fstream file;
                     std::string newFilename;
                     newFilename.append("smearedEnergyDensity_");
                     newFilename.append(_filename);
                     file.open(newFilename, std::ios::out);

                     if(!file.is_open()){return;}

                     for(int x = 0; x < grid -> getMaxSitesPerDirection(); x += accuracy){
                            for(int y = 0; y < grid -> getMaxSitesPerDirection(); y += accuracy){

                            Site site(x,y);
                            file << ((double) x / 100.0) - 15.0 << "\t" << ((double) y / 100.0) - 15.0
                            << "\t" << grid -> getSite(site) << std::endl;


                            }
                     }

                     file.close();
              }
};

template<class floatT>
class EnergyDensity{
       private:
              Grid<floatT> _grid;
              Grid<floatT> _gridSmeared;
              floatT _NNCross;
              floatT _RadNucleon;

       public:
              // constructor
              EnergyDensity(floatT newNNCross, int elems) : _grid(elems) , _gridSmeared(elems),
                                   _NNCross(newNNCross), _RadNucleon(std::sqrt(newNNCross)*0.178412*100) {}
                                                         // the factor 0.178412 is a conversion factor
                                                         // from cross section to radius in fm.
                                                         // the factor 100 is necessary due to the grid size



              // calculate energy density in MeV from NColl
              // (15.0 -> C. Schmidt, 800.0 -> Dinner N. Borghini & ???)
              void energyDensity(Grid<floatT> & grid){

                     floatT EDens = 0;

                     for (int y = 0; y < _grid.getMaxSitesPerDirection(); y++) {
                            for (int x = 0; x < _grid.getMaxSitesPerDirection(); x++) {
                                   Site site(x,y);
                                   int NColl = (int) grid.getSite(site);

                                   // 15 * 800 * NColl^4 / 46^4
              	              EDens = 0.002680093338717343 * NColl * NColl * NColl * NColl;

                                   _grid.setSite(site,EDens);
                            }
                     }
              }

              Grid<floatT> * getSmearedEnergyDensData(){
                     return & _gridSmeared;
              }

              void smearedEnergyDensity(){
                     for (int y = 0; y < _grid.getMaxSitesPerDirection(); y++) {
                            for (int x = 0; x < _grid.getMaxSitesPerDirection(); x++) {
                                   Site site(x, y);
                                   if (_grid.getSite(site) != 0 ) {
                                          for (int ynew = y - (int) _RadNucleon; ynew <= y + (int) _RadNucleon; ynew++) {
                                                 for (int xnew = x - (int) _RadNucleon; xnew <= x + (int) _RadNucleon; xnew++) {
                                                               Site siteNew(xnew,ynew);
                                                               _gridSmeared.addSite(siteNew,_grid.getSite(site));
                                                 }
                                          }
                                   }
                            }
                     }
              }
};


int main(int argc, char const *argv[]) {
       std::clock_t start;
       double duration;

       // number of events
       int NEvents;

       std::cout << "Number of events: \n";

       std::cin >> NEvents;

       // number of nucleons in the core of the element, e.g. 208 for Pb
       int NNucleonsCore;

       std::cout << "Number of nucleons per core: \n";

       std::cin >> NNucleonsCore;

       // nucleon nucleon cross section for nucleon radius
       PREC NNCross;

       std::cout << "Nucleon Nucleon cross section in mb: \n";

       std::cin >> NNCross;

       start = std::clock();

       int elems = 3001;

       std::cout << "Initializing the data grid" << '\n';

       Grid<PREC> rawDataGrid(elems);

       std::cout << "Initializing the energy density computation" << '\n';

       EnergyDensity<PREC> energDens(NNCross, rawDataGrid.getMaxSitesPerDirection());

       std::cout << "Read data " << '\n';

       std::string filename = "Pb67.6.txt";

       FileWriter<PREC> file(NEvents,NNucleonsCore, filename);

       file.readFile(rawDataGrid);

       std::cout << "Compute energy density" << '\n';

       energDens.energyDensity(rawDataGrid);

       std::cout << "Smear energy density" << '\n';

       energDens.smearedEnergyDensity();

       file.writeFileEDens(energDens.getSmearedEnergyDensData());

       duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

       std::cout<< duration <<'\n';

       return 0;
}





// Testing routins

// Test indexer:
void indexerTest(int max){
       Grid<PREC> grid(max);

       std::cout << "xyz" << '\n';

       //loop over all sites
       for (size_t y = 0; y <= max; y++) {
              for (size_t x = 0; x <= max; x++) {
                     Site si(x,y);
                     grid.setSite(si, grid.getIndex(si));
              }
       }


       //loop over all sites
       for (size_t y = 0; y <= max; y++) {
              for (size_t x = 0; x <= max; x++) {
                     Site si(x,y);
                     std::cout << x << y << "\t" << grid.getIndex(si)
                            - grid.getSite(si) << '\n';

              }
       }

}
