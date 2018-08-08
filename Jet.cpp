#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<cmath>
#include<complex>
#include <vector>

#define PREC double

struct Site;
template<class floatT> class GridAccessor;
template<class floatT> class Grid;

struct Site{
       std::vector<int> VCoordinates;

       // constructor
       Site(int x, int y, int z){
              VCoordinates.push_back(x);
              VCoordinates.push_back(y);
              VCoordinates.push_back(z);
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
                     return site.VCoordinates[0] + site.VCoordinates[1] + site.VCoordinates[2] + site.VCoordinates[1]*max
                                   + site.VCoordinates[2]*(2*max + max*max);
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
              Grid(int maxSitesPerDirection) : _VData(
                                     3*maxSitesPerDirection
                                   + 3*maxSitesPerDirection*maxSitesPerDirection
                                   + maxSitesPerDirection * maxSitesPerDirection * maxSitesPerDirection
                                   + 1 ){

                     _maxSitesPerDirection = maxSitesPerDirection;
              }

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
              GridAccessor<floatT> & _gridAcc;
       public:
              // constructor
              FileWriter(GridAccessor<floatT> & newGridAcc, int newNEvents, int newNNucleonsCore)
                                          : _gridAcc(newGridAcc) {
                     _NEvents = newNEvents;
                     _NNucleonsCore = newNNucleonsCore;
              }

              void readFile(std::string filename){
                     std::fstream file;
                     file.open(filename.c_str(), std::ios::in);

                     if(!file.is_open()){return;}

                     std::vector<floatT> row;

                     // loop through the data file which format is clarified by
                     // x corrd \t y coord \t z coord \t NColl
                     for(int i = 0; i < 2 * _NEvents * _NNucleonsCore; ++i){
                            for(int col = 0; col < 4 ; ++col){
                                   file >> row[col];
                            }
                            int x = (row[0] + 15)*100;
                            int y = (row[1] + 15)*100;
                            int z = (row[2] + 15)*100;

                            if(x > _gridAcc.getMaxSitesPerDirection() ||
                                   y > _gridAcc.getMaxSitesPerDirection() ||
                                   z > _gridAcc.getMaxSitesPerDirection()){
                                          std::cout << "ERROR@readFile: Coordinates out of grid!" << '\n';
                                   }


                            // compute the site
                            Site site(x, y, z);

                            // set NColl on the grid
                            _gridAcc.setSite(site, row[3]);
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
              GridAccessor<floatT> & _gridAcc;
              Grid<floatT> _grid;
              Grid<floatT> _gridSmeared;
              floatT _NNCross;

       public:
              // constructor
              EnergyDensity(GridAccessor<floatT> & gridAcc, floatT newNNCross): _gridAcc(gridAcc){
                     Grid<floatT> _grid(gridAcc.getMaxSitesPerDirection());
                     Grid<floatT> _gridSmeared(gridAcc.getMaxSitesPerDirection());
                     _NNCross = newNNCross;
              }

              // calculate energy density in MeV from NColl (15.0 -> C. Schmidt, 800.0 -> Dinner N. Borghini & ???)
              void energyDensity(Site site){
                     int NColl = (int) _gridAcc.getSite(site);

              	floatT EDens = 15.0*((floatT) NColl/46.0)*((floatT) NColl/46.0)
                                   *((floatT) NColl/46.0)*((floatT) NColl/46.0)*800.0;

                     _grid.getAccsessor().setSite(site,EDens);
              }

              //// nucleon radius in fm from cross section
              //floatT RadiusNucleon(floatT NNCross){
              //       floatT RadNucleon = std::sqrt(NNCross)*0.178412;
              //       return RadNucleon;
              //}
};


int main(int argc, char const *argv[]) {
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

       int elems = 4;

       Grid<PREC> grid(elems);

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
