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

// structure to define the lattice points leaving out the z component due to the contemplated event plane
struct Site{
       std::vector<int> VCoordinates;

       // constructor
       Site(int x, int y){
              VCoordinates.push_back(x);
              VCoordinates.push_back(y);
       }
};

// class to accsess the data
template<class floatT>
class GridAccessor{
       private:
              Grid<floatT> & _grid;
       public:
              GridAccessor(Grid<floatT> & newGrid) : _grid(newGrid){}

              int getIndex(Site site){
                     int max = _grid._maxSitesPerDirection;
                     return site.VCoordinates[0] + site.VCoordinates[1] + site.VCoordinates[1]*max;
              }

              void setSite(Site site, floatT Ncoll){
                     _grid._VData[getIndex(site)] = Ncoll;
              }

              floatT getSite(Site site){
                     return _grid._VData.at(getIndex(site));
              }

              int getMaxSitesPerDirection(){
                     return _grid._maxSitesPerDirection;
              }
};

// class to store the data
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

              GridAccessor<floatT> getAccsessor(){
                     GridAccessor<floatT> newGridAccessor(*this);
                     return newGridAccessor;
              }

       friend class GridAccessor<floatT>;
};


template<class floatT>
class FileWriter{
       private:
              int _NEvents;
              int _NNucleonsCore;
       public:
              // constructor
              FileWriter(int newNEvents, int newNNucleonsCore) : _NEvents(newNEvents) ,
                                                        _NNucleonsCore(newNNucleonsCore) {}

              void readFile(GridAccessor<floatT> gridAcc ,std::string filename){
                     std::fstream file;
                     file.open(filename.c_str(), std::ios::in);

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

                            if(x > gridAcc.getMaxSitesPerDirection() ||
                               y > gridAcc.getMaxSitesPerDirection()  ) {
                                          std::cout << "ERROR@readFile: Coordinates out of grid!" << '\n';
                            }

                            // compute the site
                            Site site(x, y);

                            // set NColl on the grid
                            gridAcc.setSite(site, row[3]);
                     }

                     file.close();
              }

              //void writeFileEDens(std::string filename, ){
              //       std::fstream file;
              //       file.open(filename.c_str(), std::ios::out);

              //       if (file.is_open())
              //       {
              //              file << "This is another line.\n";
              //              file.close();
              //       }
              //       else cout << "Unable to open file";
              //}
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
              void energyDensity(GridAccessor<floatT> gridAcc){

                     floatT EDens = 0;

                     for (int y = 0; y < _grid.getAccsessor().getMaxSitesPerDirection(); y++) {
                            for (int x = 0; x < _grid.getAccsessor().getMaxSitesPerDirection(); x++) {
                                   Site site(x,y);
                                   int NColl = (int) gridAcc.getSite(site);

                                   // 15 * 800 * NColl^4 / 46^4
              	              EDens = 0.002680093338717343 * NColl * NColl * NColl * NColl;

                                   _grid.getAccsessor().setSite(site,EDens);
                            }
                     }
              }

              void smearedEnergyDensity(GridAccessor<floatT> gridAcc){
                     for (int y = 0; y < _grid.getAccsessor().getMaxSitesPerDirection(); y++) {
                            for (int x = 0; x < _grid.getAccsessor().getMaxSitesPerDirection(); x++) {
                                   Site site(x, y);
                                   if (gridAcc.getSite(site) != 0 ) {
                                          for (int ynew = 0; ynew < _gridSmeared.getAccsessor().getMaxSitesPerDirection(); ynew++) {
                                                 for (int xnew = 0; xnew < _gridSmeared.getAccsessor().getMaxSitesPerDirection(); xnew++) {
                                                        if (_RadNucleon > std::hypot((x - xnew),(y - ynew))) {
                                                               Site siteNew(xnew,ynew);
                                                               _gridSmeared.getAccsessor().setSite(siteNew,_grid.getAccsessor().getSite(site));
                                                        }
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

       start = std::clock();

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

       int elems = 3001;

       std::cout << "Initializing the data grid" << '\n';

       Grid<PREC> rawDataGrid(elems);

       std::cout << "Initializing the energy density computation" << '\n';

       EnergyDensity<PREC> energDens(NNCross, rawDataGrid.getAccsessor().getMaxSitesPerDirection());

       std::cout << "Read data " << '\n';

       FileWriter<PREC> file(NEvents,NNucleonsCore);

       std::string filename = "Test.dat";

       file.readFile(rawDataGrid.getAccsessor(), filename);

       std::cout << "Compute energy density" << '\n';

       energDens.energyDensity(rawDataGrid.getAccsessor());

       duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

       std::cout<< duration <<'\n';

       return 0;
}




/*
// Testing routins

// Test indexer:
void indexerTest(int max){
       Grid<PREC> grid(max);

       std::cout << "xyz" << '\n';

       //loop over all sites
       for (size_t z = 0; z <= max; z++) {
              for (size_t y = 0; y <= max; y++) {
                     for (size_t x = 0; x <= max; x++) {
                            Site si(x,y,z);
                            grid.getAccsessor().setSite(si, grid.getAccsessor().getIndex(si));
                     }
              }
       }

       //loop over all sites
       for (size_t z = 0; z <= max; z++) {
              for (size_t y = 0; y <= max; y++) {
                     for (size_t x = 0; x <= max; x++) {
                            Site si(x,y,z);
                            std::cout << x << y << z << "\t" << grid.getAccsessor().getIndex(si)
                                   - grid.getAccsessor().getSite(si) << '\n';

                     }
              }
       }
}


*/
