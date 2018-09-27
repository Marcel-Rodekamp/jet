#include<iostream>
#include<iomanip>
#include<fstream>
#include<stdlib.h>
#include<cmath>
#include<complex>
#include<vector>
#include<omp.h>

// floats are not available (not working with class complex)
#define PI 3.14159265358979323846
#define Steps 500
#define PREC double

class Site;
template<class floatT> class Grid;
template<class floatT> class FileWriter;
template<class floatT> class EnergyDensity;
template<class floatT> class Eccentricity;

void testIndexer(int max);
void testRunThroughGrid(int max);
void testReadWrite();
void testFileReader();

// structure to define the lattice points leaving out the z component due to the contemplated event plane
class Site{
       private:
              std::vector<int> VCoordinates;

       public:
              // constructor
              Site(int x, int y): VCoordinates(2){
                     VCoordinates.at(0) = x;
                     VCoordinates.at(1) = y;
              }

              int inline x(){
                     return VCoordinates.at(0);
              }

              int inline y(){
                     return VCoordinates.at(1);
              }

              void inline xUp(){
                     VCoordinates.at(0) += 1;
              }

              void inline yUp(){
                     VCoordinates.at(1) += 1;
              }

              void inline xDown(){
                     VCoordinates.at(0) -= 1;
              }

              void inline yDown(){
                     VCoordinates.at(1) -= 1;
              }

              void inline setX(int newX){
                     VCoordinates.at(0) = newX;
              }

              void inline setY(int newY){
                     VCoordinates.at(1) = newY;
              }

              void inline printSite(){
                     std::cout << "(" << this -> x() << "," << this -> y() << ") \t";
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
              // constructor
              Grid(int maxSitesPerDirection) :
                                   _VData( 2*maxSitesPerDirection
                                          + maxSitesPerDirection*maxSitesPerDirection
                                          + 1 ) ,
                                   _maxSitesPerDirection(maxSitesPerDirection) {}

              int inline getIndex(Site site){
                     return site.x() + site.y() + site.y() * _maxSitesPerDirection;
              }

              void inline setSite(Site site, floatT value){
                     _VData.at(getIndex(site)) = value;
              }

              floatT inline getSite(Site site){
                     return _VData.at(getIndex(site));
              }

              void inline addSite(Site site, floatT value){
                     _VData.at(getIndex(site)) += value;
              }

              int inline getMaxSitesPerDirection(){
                     return _maxSitesPerDirection;
              }

              // run through the whole grid
              bool runThroughGrid(Site & site){
                     // compute possible new x value
                     int newCoordX = site.x() + 1;

                     // check if the x value is on the grid
                     if(newCoordX <= _maxSitesPerDirection){
                            site.xUp();
                     }
                     else{
                            // if the x value is not on the grid compute possible new y value
                            int newCoordY = site.y() + 1;

                            // check if the y value is on  the grid
                            if(newCoordY <= _maxSitesPerDirection){
                                   site.yUp();
                                   site.setX(0);
                            }
                            else{
                                   // end the loop over the grid
                                   return false;
                            }
                     }

                     return true;
              }
};

// class to write data in files
template<class floatT>
class FileWriter{
       private:
                  const std::string _fileName;
                  std::ofstream _fileStream;

                  // initializes class, automatically called by constructors
                  void init() {
                         _fileStream.open(_fileName.c_str());

                          if(!_fileStream.is_open()){
                                 std::cerr << "Error@FileWriter:File could not be opened" << '\n';
                                 return;
                          }
                  }

       public:

              // constructor
              FileWriter(std::string fname) :
                            _fileName(fname) {
                     init();
              }

              std::ofstream & operator<< (const floatT &obj) {
                     _fileStream << obj;
                     return _fileStream;
              }

              std::ofstream & operator<< (const std::string &string_obj) {
                     _fileStream << string_obj;
                     return _fileStream;
              }

              // Close the ostream
              ~FileWriter() {
                     _fileStream.close();
              }
};

// class to read data from file
template<class floatT>
class FileReader{
       private:
           const std::string _fileName;
           std::ifstream _fileStream;
           std::vector<int> _data;

           // initializes class, automatically called by constructors
           void init() {
                  _fileStream.open(_fileName.c_str(),std::ios::in);

                   if(!_fileStream.is_open()){
                          std::cerr << "Error@FileReader:File could not be opened" << '\n';
                          return;
                   }
           }

    public:
       // constructor
       FileReader(int NEvents, int NNucleonsCore, std::string fname) :
                     _fileName(fname),
                     // 2 * Nucleons * Events * Values to store
                     _data(3 * 2 * NNucleonsCore * NEvents){

              //init();
       };

       void readData(int NEvents, int NNucleonsCore){
               std::ifstream _fileStream;
                _fileStream.open(_fileName.c_str(),std::ios::in);

              std::vector<floatT> row(4);

              // loop through the data file which format is clarified by
              // x corrd \t y coord \t z coord \t NColl
              for(int i = 0; i < 2 * 3 * NEvents * NNucleonsCore; i += 3){
                     for(int col = 0; col < 4 ; col++) {
                            _fileStream >> row.at(col);
                     }

                     // translate the euclidean coordinates into the grid coordinates
                     // leaving out the z component
                     int x = (row.at(0) + 15)*100;
                     int y = (row.at(1) + 15)*100;

                     _data.at(i) = x;
                     _data.at(i+1) = y;
                     _data.at(i+2) = row.at(3);

                     row.at(0) = 0.;
                     row.at(1) = 0.;
                     row.at(2) = 0.;
                     row.at(3) = 0.;
              }
       }
       void getData(int event, int NNucleonsCore, Grid<floatT> & grid){
              Site site(0,0);
              for(int i = 3 * 2 * (event - 1) * NNucleonsCore; i < 3 * 2 * event * NNucleonsCore; i += 3){
                     site.setX(_data.at(i));
                     site.setY(_data.at(i+1));

                     if(_data.at(i) > grid.getMaxSitesPerDirection() || _data.at(i+1) > grid.getMaxSitesPerDirection()){
                            std::cerr << "Error@getData: Nucleon out of grid with: " << _data.at(i+2) << '\n';
                            continue;
                     }

                     grid.setSite(site, _data.at(i+2));

              }
       }

       ~FileReader(){
              _fileStream.close();
       }

};

// class to calculate the energy density at a grid point and smear it out
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
                                   _NNCross(newNNCross), _RadNucleon((std::sqrt(newNNCross)*0.178412/2.)*100) {}
                                                         // factor 1/2 convention from [arXiv:1710.07098]
                                                         // the factor 0.178412 is a conversion factor
                                                         // from cross section to radius in fm.
                                                         // the factor 100 is necessary due to the grid steps

              // calculate energy density in MeV from NColl
              // (15.0 -> C. Schmidt, 800.0 in hot spots -> N. Borghini & U. Heinz)
              void energyDensity(Grid<floatT> & grid){
                     floatT EDens = 0;
                     Site site(0,0);
                     int NColl;

                     do{
                            NColl = (int) grid.getSite(site);

                            // (15*800/46^4) * NColl^4
                            EDens = 0.002680093338717343 * NColl * NColl * NColl * NColl;

                            _grid.setSite(site,EDens);
                     }while(grid.runThroughGrid(site));
              }

              void smearedEnergyDensity(){
                     Site site(0, 0);
                     Site tmpSite(0,0);
                     // TODO this is not so nice!
                     const int tmpRadNucleon = (int) _RadNucleon;

                     do{
                            if (_grid.getSite(site) != 0 ) {
                                   // loop through a square with side length tmpRadNucleon centered on the
                                   // not zero site
                                   // +3 and -3 to prevent edge effects
                                   for(int x = site.x() - tmpRadNucleon - 3; x < site.x() + tmpRadNucleon + 3; x++){
                                          if(x < 0){
                                                 std::cerr << "Error@smearedEnergyDensity: x value to small" << '\n';
                                                 return ;
                                          }
                                          if(x > _gridSmeared.getMaxSitesPerDirection()){
                                                 std::cerr << "Error@smearedEnergyDensity: x value to big" << '\n';
                                                 return ;
                                          }
                                   for(int y = site.y() - tmpRadNucleon - 3; y < site.y() + tmpRadNucleon + 3; y++){
                                          if(y < 0){
                                                 std::cerr << "Error@smearedEnergyDensity: y value to small" << '\n';
                                                 return ;
                                          }
                                          if(y > _gridSmeared.getMaxSitesPerDirection()){
                                                 std::cerr << "Error@smearedEnergyDensity: y value to big" << '\n';
                                                 return ;
                                          }
                                          if(std::abs((x-site.x())*(x-site.x()) + (y-site.y())*(y-site.y())) < _RadNucleon*_RadNucleon){
                                                 tmpSite.setX(x);
                                                 tmpSite.setY(y);
                                                 _gridSmeared.addSite(tmpSite,_grid.getSite(site));
                                          }
                                   }
                                   }
                            }
                     }while( _grid.runThroughGrid(site));

              }

              Grid<floatT> * getSmearedEnergyDensData(){
                     return & _gridSmeared;
              }

              friend class IntegratedEnergyDensity<floatT>;
              friend class MeanStartingPoint<floatT>;
              friend class Eccentricity<floatT>;
};

// class to compute the eccentricities e2 - e6
template<class floatT>
class Eccentricity{
       private:
              std::vector<floatT> _EventEccentricity;

              Grid<floatT> & _smearEnerDensGrid;

              std::vector<std::complex<floatT>> _CounterAngle;
              std::vector<floatT> _DenomAngle;

              floatT _Ecc2;
              floatT _Ecc3;
              floatT _Ecc4;
              floatT _Ecc5;
              floatT _Ecc6;

              floatT _meanX;
              floatT _meanY;

              // function for the calculation of the mean x value
              // trapezoidal integration 2D of (int dxdy x*e(x,y))/(int dxdy e(x,y))
              void inline _calculateMeanX(){
                     floatT IntNominator = 0.;
                     floatT IntDenominator = 0.;
                     for (int x = 0; x < _smearEnerDensGrid.getMaxSitesPerDirection(); x++) {
                            for (int y = 0; y < _smearEnerDensGrid.getMaxSitesPerDirection(); y++) {
                                   Site tmpSite(x,y);
                                   if (_smearEnerDensGrid.getSite(tmpSite) != 0.0) {
                                          // conversion to real coordinate system
                                          floatT realX = tmpSite.x()/100. - 15.;
                                          IntNominator += (realX * _smearEnerDensGrid.getSite(tmpSite));
                                          IntDenominator += _smearEnerDensGrid.getSite(tmpSite);
                                   }
                            }
                     }
                     _meanX = IntNominator / IntDenominator;
              }

              // function for the calculation of the mean y value
              // trapezoidal integration 2D of (int dxdy y*e(x,y))/(int dxdy e(x,y))
              void inline _calculateMeanY(){
                     floatT IntNominator = 0.;
                     floatT IntDenominator = 0.;
                     for (int x = 0; x < _smearEnerDensGrid.getMaxSitesPerDirection(); x++) {
                            for (int y = 0; y < _smearEnerDensGrid.getMaxSitesPerDirection(); y++) {
                                   Site tmpSite(x,y);
                                   if (_smearEnerDensGrid.getSite(tmpSite) != 0.0) {
                                          // conversion to real coordinate system
                                          floatT realY = tmpSite.y()/100. - 15.;
                                          IntNominator += (realY * _smearEnerDensGrid.getSite(tmpSite));
                                          IntDenominator += _smearEnerDensGrid.getSite(tmpSite);
                                   }
                            }
                     }
                     _meanY = IntNominator / IntDenominator;
              }

              // compute the ecc2
              void inline _computeEcc2() {
                     // define the imaginary I
                     const std::complex<floatT> I(0.0, 1.0);
                     std::complex<floatT> IntNominator = (0.0);
                     floatT IntDenominator = 0.;

                     for (int x = 0; x < _smearEnerDensGrid.getMaxSitesPerDirection(); x++) {
                            for (int y = 0; y < _smearEnerDensGrid.getMaxSitesPerDirection(); y++) {
                                   Site tmpSite2(x,y);
                                   // conversion to real coordinate system shifted by mean values
                                   floatT realX = (tmpSite2.x()/100. - 15.) - _meanX;
                                   floatT realY = (tmpSite2.y()/100. - 15.) - _meanY;

                                   // calculated site which is shifted by the mean x,y
                                   int xs = (int) ((realX + _meanX + 15.)*100.);
                                   int ys = (int) ((realY + _meanY + 15.)*100.);
                                   Site tmpSite1(xs,ys);

                                   if (_smearEnerDensGrid.getSite(tmpSite1) != 0.0) {
                                          IntNominator += (realX*realX - realY*realY + 2.0*I*realX*realY)
                                                        * _smearEnerDensGrid.getSite(tmpSite1);

                                          IntDenominator += (realX*realX + realY*realY)
                                                        * _smearEnerDensGrid.getSite(tmpSite1);
                                   }
                            }
                     }
                     _Ecc2 = std::abs(IntNumerator / IntDenominator);
              }

              // compute the ecc3
              void inline _computeEcc3() {
                     // define the imaginary I
                     const std::complex<floatT> I(0.0, 1.0);
                     std::complex<floatT> IntNumerator = (0.0);
                     floatT IntDenominator = 0.;

                     for (int x = 0; x < _smearEnerDensGrid.getMaxSitesPerDirection(); x++) {
                            for (int y = 0; y < _smearEnerDensGrid.getMaxSitesPerDirection(); y++) {
                                   Site tmpSite2(x,y);
                                   // conversion to real coordinate system shifted by mean values
                                   floatT realX = (tmpSite2.x()/100. - 15.) - _meanX;
                                   floatT realY = (tmpSite2.y()/100. - 15.) - _meanY;

                                   // calculated site which is shifted by the mean x,y
                                   int xs = (int) ((realX + _meanX + 15.)*100.);
                                   int ys = (int) ((realY + _meanY + 15.)*100.);
                                   Site tmpSite1(xs,ys);

                                   if (_smearEnerDensGrid.getSite(tmpSite1) != 0.0) {
                                          IntNumerator += (realX*realX*realX - 3.0*realX*realY*realY
                                                 + I*(3.0*realX*realX*realY - realY*realY*realY))
                                                 * _smearEnerDensGrid.getSite(tmpSite1);

                                          IntDenominator += ((realX*realX + realY*realY)
                                                        * std::sqrt(realX*realX + realY*realY))
                                                        * _smearEnerDensGrid.getSite(tmpSite1);
                                   }
                            }
                     }
                     _Ecc3 = std::abs(IntNumerator / IntDenominator);
              }

              // compute the ecc4
              void inline _computeEcc4() {
                     // define the imaginary I
                     const std::complex<floatT> I(0.0, 1.0);
                     std::complex<floatT> IntNumerator = (0.0);
                     floatT IntDenominator = 0.;

                     for (int x = 0; x < _smearEnerDensGrid.getMaxSitesPerDirection(); x++) {
                            for (int y = 0; y < _smearEnerDensGrid.getMaxSitesPerDirection(); y++) {
                                   Site tmpSite2(x,y);
                                   // conversion to real coordinate system shifted by mean values
                                   floatT realX = (tmpSite2.x()/100. - 15.) - _meanX;
                                   floatT realY = (tmpSite2.y()/100. - 15.) - _meanY;

                                   // calculated site which is shifted by the mean x,y
                                   int xs = (int) ((realX + _meanX + 15.)*100.);
                                   int ys = (int) ((realY + _meanY + 15.)*100.);
                                   Site tmpSite1(xs,ys);

                                   if (_smearEnerDensGrid.getSite(tmpSite1) != 0.0) {
                                          IntNumerator += (realX*realX*realX*realX
                                                        + realY*realY*realY*realY
                                                        - 6.0*realX*realX*realY*realY
                                                        + 4.0*I*(realX*realX*realX*realY
                                                        - realX*realY*realY*realY))
                                                        * _smearEnerDensGrid.getSite(tmpSite1);

                                          IntDenominator += (realX*realX*realX*realX
                                                               + realY*realY*realY*realY
                                                               + 2.0*realX*realX*realY*realY)
                                                               * _smearEnerDensGrid.getSite(tmpSite1);
                                   }
                            }
                     }
                     _Ecc4 = std::abs(IntNumerator / IntDenominator);
              }

              // compute the ecc5
              void inline _computeEcc5() {
                     // define the imaginary I
                     const std::complex<floatT> I(0.0, 1.0);
                     std::complex<floatT> IntNumerator = (0.0);
                     floatT IntDenominator = 0.;

                     for (int x = 0; x < _smearEnerDensGrid.getMaxSitesPerDirection(); x++) {
                            for (int y = 0; y < _smearEnerDensGrid.getMaxSitesPerDirection(); y++) {
                                   Site tmpSite2(x,y);
                                   // conversion to real coordinate system shifted by mean values
                                   floatT realX = (tmpSite2.x()/100. - 15.) - _meanX;
                                   floatT realY = (tmpSite2.y()/100. - 15.) - _meanY;

                                   // calculated site which is shifted by the mean x,y
                                   int xs = (int) ((realX + _meanX + 15.)*100.);
                                   int ys = (int) ((realY + _meanY + 15.)*100.);
                                   Site tmpSite1(xs,ys);

                                   if (_smearEnerDensGrid.getSite(tmpSite1) != 0.0) {
                                          IntNumerator += (realX*realX*realX*realX*realX
                                                 - 10.0*realX*realX*realX*realY*realY
                                                 + 5.0*realX*realY*realY*realY*realY
                                                 + I*(realY*realY*realY*realY*realY
                                                 + 5.0*realX*realX*realX*realX*realY
                                                 - 10.0*realX*realX*realY*realY*realY))
                                                 * _smearEnerDensGrid.getSite(tmpSite1);

                                          IntDenominator += ((realX*realX + realY*realY)
                                                               * (realX*realX + realY*realY)
                                                               * std::sqrt(realX*realX + realY*realY))
                                                               * _smearEnerDensGrid.getSite(tmpSite1);
                                   }
                            }
                     }
                     _Ecc5 = std::abs(IntNumerator / IntDenominator);
              }

              // compute the ecc6
              void inline _computeEcc6() {
                     // define the imaginary I
                     const std::complex<floatT> I(0.0, 1.0);
                     std::complex<floatT> IntNumerator = (0.0);
                     floatT IntDenominator = 0.;

                     for (int x = 0; x < _smearEnerDensGrid.getMaxSitesPerDirection(); x++) {
                            for (int y = 0; y < _smearEnerDensGrid.getMaxSitesPerDirection(); y++) {
                                   Site tmpSite2(x,y);
                                   // conversion to real coordinate system shifted by mean values
                                   floatT realX = (tmpSite2.x()/100. - 15.) - _meanX;
                                   floatT realY = (tmpSite2.y()/100. - 15.) - _meanY;

                                   // calculated site which is shifted by the mean x,y
                                   int xs = (int) ((realX + _meanX + 15.)*100.);
                                   int ys = (int) ((realY + _meanY + 15.)*100.);
                                   Site tmpSite1(xs,ys);

                                   if (_smearEnerDensGrid.getSite(tmpSite1) != 0.0) {
                                          IntNumerator += (realX*realX*realX*realX*realX*realX
                                                 - realY*realY*realY*realY*realY*realY
                                                 - 15.0*realX*realX*realX*realX*realY*realY
                                                 + 15.0*realX*realX*realY*realY*realY*realY
                                                 + I*(6.0*realX*realX*realX*realX*realX*realY
                                                 - 20.0*realX*realX*realX*realY*realY*realY
                                                 + 6.0*realX*realY*realY*realY*realY*realY))
                                                 * _smearEnerDensGrid.getSite(tmpSite1);

                                          IntDenominator += (realX*realX*realX*realX*realX*realX
                                                        + realY*realY*realY*realY*realY*realY
                                                        + 3.0*realX*realX*realY*realY*realY*realY
                                                        + 3.0*realX*realX*realX*realX*realY*realY)
                                                        * _smearEnerDensGrid.getSite(tmpSite1);
                                   }
                            }
                     }
                     _Ecc6 = std::abs(IntNumerator / IntDenominator);
              }

       public:
              // constructor
              Eccentricity(EnergyDensity<floatT> & newEnergDens):
                            _EventEccentricity(5),
                            _smearEnerDensGrid(* newEnergDens.getSmearedEnergyDensData()),
                            _CounterAngle(5),
                            _DenomAngle(5),
                            _Ecc2(),
                            _Ecc3(),
                            _Ecc4(),
                            _Ecc5(),
                            _Ecc6(),
                            _meanX(),
                            _meanY(){}

              // calculate the mean x and y values of e(x,y) and the eccentricities e2 ... e6 for one event
              void computeEccentricity(){
                     _calculateMeanX();
                     _calculateMeanY();
                     _computeEcc2();
                     _computeEcc3();
                     _computeEcc4();
                     _computeEcc5();
                     _computeEcc6();
              }

              // functions to get the values of the eccentricities for one event
              floatT getEcc2() {
                     return _Ecc2;
              }

              floatT getEcc3() {
                     return _Ecc3;
              }

              floatT getEcc4() {
                     return _Ecc4;
              }

              floatT getEcc5() {
                     return _Ecc5;
              }

              floatT getEcc6() {
                     return _Ecc6;
              }
};

int main(int argc, char const *argv[]) {

       // read number of events
       int NEvents;
       std::cout << "Number of events: \n";
       std::cin >> NEvents;

       // read number of nucleons in the core of the element, e.g. 208 for Pb
       int NNucleonsCore;
       std::cout << "Number of nucleons per core: \n";
       std::cin >> NNucleonsCore;

       // read nucleon nucleon cross section for nucleon radius
       PREC NNCross;
       std::cout << "Nucleon Nucleon cross section in mb: \n";
       std::cin >> NNCross;

       // define size of grid = elems*elems + 2 elems for "axes"
       int elems = 3001;

       // build the FileReader class
       FileReader<PREC> fr(NEvents, NNucleonsCore, "Pb67.6.txt");

       // read data from file
       fr.readData(NEvents, NNucleonsCore);

       // file for the eccentricities
       std::string filename2 = "Eccentricities.dat";
       FileWriter<PREC> file2(filename2);

       #pragma omp parallel for
       for (int Event = 1; Event <= NEvents; Event++) {
              std::cout << Event << '\n';

              // build grid class
              Grid<PREC> rawDataGrid(elems);

              // build energy density class
              EnergyDensity<PREC> energDens(NNCross, elems);

              #pragma omp critical
              // get data from FileReader class
              fr.getData(Event, NNucleonsCore, rawDataGrid);

              // calculate the energy density from the data grid
              energDens.energyDensity(rawDataGrid);

              // smear the energy density over a predefined radius (in EnergyDensity class)
              energDens.smearedEnergyDensity();

              // build the Eccentricity class
              Eccentricity<PREC> ecc(energDens);

              // compute the eccentricities
              ecc.computeEccentricity();

              // write eccentricities in "Eccentricities.dat"
              #pragma omp critical
              file2  << ecc.getEcc2() << "\t"
                     << ecc.getEcc3() << "\t"
                     << ecc.getEcc4() << "\t"
                     << ecc.getEcc5() << "\t"
                     << ecc.getEcc6() << std::endl;

              // write one example smeared energy density, only 1/10 of the values --> smaller data file
              if(Event == 1){
                     std::string filename1 = "EDensity.dat";
                     FileWriter<PREC> file1(filename1);
                     // write smeared Energy density in "EDensity.dat"
                     for (int x = 0; x < energDens.getSmearedEnergyDensData() -> getMaxSitesPerDirection(); x += 10) {
                            for (int y = 0; y < energDens.getSmearedEnergyDensData() -> getMaxSitesPerDirection(); y += 10) {
                                   Site site(x,y);
                                   file1  << site.x()/100. - 15. << "\t"
                                          << site.y()/100. - 15. << "\t"
                                          << energDens.getSmearedEnergyDensData() -> getSite(site)
                                          << std::endl;
                            }
                     }
              }
       }
       return 0;
}

// testing routines
// test indexer:
void testIndexer(int max){
       Grid<PREC> grid(max);

       std::cout << "xyz" << '\n';

       //loop over all sites
       for (int y = 0; y <= max; y++) {
              for (int x = 0; x <= max; x++) {
                     Site si(x,y);
                     grid.setSite(si, grid.getIndex(si));
              }
       }


       //loop over all sites
       for (int y = 0; y <= max; y++) {
              for (int x = 0; x <= max; x++) {
                     Site si(x,y);
                     std::cout << x << y << "\t" << grid.getIndex(si)
                            - grid.getSite(si) << '\n';

              }
       }

}

// test runThroughGrid()
void testRunThroughGrid(int max){
       Grid<PREC> grid(max);

       std::cout << "runThroughGrid(Site)" << '\n';

       Site site(0,0);

       do{
              std::cout << "now at (x,y): (" << site.x() << "," << site.y() << ")" << '\n';
       }while(grid.runThroughGrid(site));

       std::cout << "\nloop through lattice\n" << '\n';

       for (int y = 0; y <= max; y++) {
              for (int x = 0; x <= max; x++){
                     Site si(x,y);
                     std::cout << "now at (x,y): (" << si.x() << "," << si.y() << ")" << '\n';
              }
       }
}

// test file read and write
void testReadWrite() {
       int NEvents = 200;
       int NNucleonsCore = 208;

       int elems = 3001;

       Site startSite(0,0);

       FileReader<PREC> fr(NEvents, NNucleonsCore, "Pb67.6.txt");
       fr.readData(NEvents, NNucleonsCore);

       std::string dataReadWrite = "ReadWrite.dat";
       FileWriter<PREC> writeData(dataReadWrite);

       for (int Event = 1; Event <= NEvents; Event++) {
              std::cout << Event << '\n';

              Grid<PREC> rawDataGrid(elems);
              fr.getData(Event, NNucleonsCore, rawDataGrid);

              for (int x = 0; x < elems; x++) {
                     for (int y = 0; y < elems; y++) {
                            Site site(x,y);
                            if(rawDataGrid.getSite(site) != 0){
                            writeData << (x/100.) - 15 << "\t" <<  (y/100.) - 15 << "\t" << rawDataGrid.getSite(site) << "\n";
                            }
                     }
              }
       }
}


// Test file reader test
void testFileReader(){
       // number of events
       int NEvents = 5;

       // number of nucleons in the core of the element, e.g. 208 for Pb
       int NNucleonsCore = 208;

       // define size of grid = elems*elems + 2 elems
       int elems = 3001;

       std::string filename = "test.dat";
       FileWriter<PREC> file(filename);
       filename = "test2.dat";
       FileWriter<PREC> file2(filename);

       Grid<PREC> rawDataGrid(elems);
       Grid<PREC> rawDataGrid2(elems);

       FileReader<PREC> data(NEvents, NNucleonsCore, "Pb67.6.txt");
       FileReader<PREC> data2(NEvents, NNucleonsCore, "Pb67.6.txt");

       data.readData(NEvents, NNucleonsCore);
       data2.readData(NEvents, NNucleonsCore);
       data2.readData(NEvents, NNucleonsCore);

       Site s3(0,0);
       do{
              if(rawDataGrid.getSite(s3) != rawDataGrid2.getSite(s3)){
                     std::cout << "(" << s3.x() << "," << s3.y() << ")" << ":\t"
                               << rawDataGrid.getSite(s3) << " "
                               << rawDataGrid2.getSite(s3) << '\n';
              }
       } while(rawDataGrid2.runThroughGrid(s3));

}
