#include<iostream>
#include<iomanip>
#include<fstream>
#include<stdlib.h>
#include<cmath>
#include<complex>
#include<vector>
#include<ctime>
#include<omp.h>

// floats are not available (not working with class complex)
#define PI 3.14159265358979323846
#define Steps 500
#define PREC double

class Site;
template<class floatT> class Grid;
template<class floatT> class FileWriter;
template<class floatT> class EnergyDensity;
template<class floatT> class IntegratedEnergyDensity;
template<class floatT> class Eccentricity;
template<class floatT> class FlowCoefficients;


void testIndexer(int max);
void testRunThroughGrid(int max);
void testCompute();
void testFileWriter(int elems, PREC NNCross, int NNucleonsCore, int NEvents);
void testFileReader();
void testReadWrite();
void testIntegratedEnergyDensity();

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
              //constructor
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

template<class floatT>
class FileWriter{
       private:
                  const std::string _fileName;
                  std::ofstream _fileStream;

                  //! initializes class, automatically called by constructors
                  void init() {
                         _fileStream.open(_fileName.c_str());

                          if(!_fileStream.is_open()){
                                 std::cerr << "Error@FileWriter:File could not be opened" << '\n';
                                 return;
                          }

                         // set high precision
                         // _fileStream.precision(15);
                         // _fileStream.setf(std::ios::scientific);
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

              //! Close the ostream
              ~FileWriter() {
                     _fileStream.close();
              }
};

template<class floatT>
class FileReader{
       private:
           const std::string _fileName;
           std::ifstream _fileStream;
           std::vector<int> _data;

           //! initializes class, automatically called by constructors
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
                                   //site.printSite();
                                   // loop through a square with side length tmpRadNucleon centered on the
                                   // not zero site
                                   for(int x = site.x() - tmpRadNucleon; x < site.x() + tmpRadNucleon; x++){
                                          if(x < 0){
                                                 std::cerr << "Error@smearedEnergyDensity: x value to small" << '\n';
                                                 return ;
                                          }
                                          if(x > _gridSmeared.getMaxSitesPerDirection()){
                                                 std::cerr << "Error@smearedEnergyDensity: x value to big" << '\n';
                                                 return ;
                                          }
                                   for(int y = site.y() - tmpRadNucleon; y < site.y() + tmpRadNucleon; y++){
                                          if(y < 0){
                                                 std::cerr << "Error@smearedEnergyDensity: y value to small" << '\n';
                                                 return ;
                                          }
                                          if(y > _gridSmeared.getMaxSitesPerDirection()){
                                                 std::cerr << "Error@smearedEnergyDensity: y value to big" << '\n';
                                                 return ;
                                          }

                                          tmpSite.setX(x);
                                          tmpSite.setY(y);
                                          _gridSmeared.addSite(tmpSite,_grid.getSite(site));
                                   }
                                   }
                            }
                     }while( _grid.runThroughGrid(site));

              }

              Grid<floatT> * getSmearedEnergyDensData(){
                     return & _gridSmeared;
              }

              friend class IntegratedEnergyDensity<floatT>;

};

//TODO _realX usw in eine Site verwandeln
//TODO get rid of the sector computation and move it on the grid -> huge speed up
template<class floatT>
class IntegratedEnergyDensity{
       private:
              // convert start point to grid coordinates
              floatT _RealX;
              floatT _RealY;
              int _JetStartX;
              int _JetStartY;

              // angle steps for radial integration
              floatT _AngleStep;

              Grid<floatT> & _smearEnerDensGrid;

              // vectors for the calculated angle and integrated energy density
              std::vector<floatT> _AngleSec1;
              std::vector<floatT> _EDensSec1;

              std::vector<floatT> _AngleSec2;
              std::vector<floatT> _EDensSec2;

              std::vector<floatT> _AngleSec3;
              std::vector<floatT> _EDensSec3;

              std::vector<floatT> _AngleSec4;
              std::vector<floatT> _EDensSec4;

              std::vector<floatT> _Integral;

              std::vector<floatT> _RadiusOne;
              std::vector<floatT> _EDensOne;
              std::vector<floatT> _RadiusTwo;
              std::vector<floatT> _EDensTwo;
              std::vector<floatT> _RadiusThree;
              std::vector<floatT> _EDensThree;
              std::vector<floatT> _RadiusFour;
              std::vector<floatT> _EDensFour;

              // calculation for the right sector (sector 1) from 7/4 pi to 1/4 pi
              void _sector1() {
                     // PhiOne is the angle with respect to the origin in real space
                     floatT PhiOne = 0.0;

                     // alpha is the angle seen from the start point of the integration
                     floatT angleAlpha = -PI / 4.0;

                     // interpolation coord
                     floatT interCoordx_1, interCoordy_1;
                     floatT interCoordx_2, interCoordy_2;
                     floatT interCoordy;

                     floatT distFromJetStart;

                     floatT Integral;

                     for(int step = 0; step < Steps; step++) {
                            // set integral to 0 for += notation
                            Integral = 0.0;

                            // calculate the slope from alpha
                     	floatT slope = std::tan(angleAlpha);

                            // point-slope form of equation for straight line y = m(x-x1)+y1
                            // calculate for each x value of the grid a value for y
                     	for(int interCoordx = _JetStartX;
                                    interCoordx < _smearEnerDensGrid.getMaxSitesPerDirection() - _RealX;
                                    interCoordx++) {

                                   interCoordy = slope * (interCoordx - _JetStartX) +  _JetStartY;

                                   // check wheather the point lies on the grid or not (needed for
                                   // integration which does not start at the origin (real (0.0,0.0))
                     		if(interCoordy < _smearEnerDensGrid.getMaxSitesPerDirection()
                                          && interCoordy >= 0.0) {


                                          // find the four gridpoints {(x1,y1),(x2,y1),(x1,y2),(x2,y2)}
                                          // for the interpolation around a point (x,y)
                     			interCoordx_1 = (int) std::floor((floatT) interCoordx);
                     			interCoordx_2 = (int) std::floor((floatT) interCoordx + 1.0);
                     			interCoordy_1 = (int) std::floor(interCoordy);
                     			interCoordy_2 = (int) std::floor(interCoordy + 1.0);

                                          // define the points where the energy density grid is read out
                     			Site siteP1(interCoordx_1, interCoordy_1);
                                          Site siteP2(interCoordx_1, interCoordy_2);
                                          Site siteP3(interCoordx_2, interCoordy_1);
                                          Site siteP4(interCoordx_2, interCoordy_2);
                                          Site siteP(interCoordx, interCoordy);

                     			floatT Interpol = _energyDensityInterpolation(
                                                 siteP1,
                                                 siteP2,
                                                 siteP3,
                                                 siteP4,
                                                 siteP);

                                          // calculate the distance from the jet origin
                     			distFromJetStart = std::hypot(
                                                 (((floatT) interCoordx - (floatT) _JetStartX) / 100.0),
                                                 ((interCoordy - (floatT) _JetStartY) / 100.0));

                                          // store point at line and the interpolated value at that point
                     			_RadiusOne.at(interCoordx - _JetStartX) = distFromJetStart;



                     			_EDensOne.at(interCoordx - _JetStartX) = Interpol;
                     		}
                     		else {
                                          // calculate the distance from the jet origin
                     			distFromJetStart = std::hypot(
                                                 (((floatT) interCoordx - (floatT) _JetStartX) / 100.0),
                                                 ((interCoordy - (floatT) _JetStartY) / 100.0));

                                          // store point at line and the interpolated value at that point
                     			_RadiusOne.at(interCoordx - _JetStartX) = distFromJetStart;
                     			_EDensOne.at(interCoordx - _JetStartX) = 0.;
                     		}
                     	}

                            // perform a trapezoidal integration along the line
                     	for(int inteStep = 0;
                                   inteStep < _smearEnerDensGrid.getMaxSitesPerDirection() - _JetStartX;
                                   inteStep++){

                                   Integral += 0.5 * (_RadiusOne.at(inteStep+1) - _RadiusOne.at(inteStep))
                                                   * (_EDensOne.at(inteStep) + _EDensOne.at(inteStep+1));
                     	}

                     // calculate the angle from the center of mass (real (0.0,0.0))
                     interCoordy = slope
                            * (_smearEnerDensGrid.getMaxSitesPerDirection() - _JetStartX - 1)
                            + (floatT) _JetStartY;

                     PhiOne = std::atan(
                            (interCoordy - ((floatT) _smearEnerDensGrid.getMaxSitesPerDirection() - 1.)/2.)
                                   /(((floatT) _smearEnerDensGrid.getMaxSitesPerDirection() - 1.)/2.));
                     if(PhiOne < 0.) {PhiOne += 2.0 * PI;}

                     // store the values for the angle and the integral vectors
                     _AngleSec1.at(step) = PhiOne;

                     _EDensSec1.at(step) = Integral;

                     angleAlpha += _AngleStep;
                     }
              }

              // calculation for the left sector (sector 2) from 3/4 pi to 5/4 pi
              void _sector2() {
                     // PhiTwo is the angle with respect to the origin in real space
                     floatT PhiTwo = 0.0;

                     // alpha is the angle seen from the start point of the integration
                     floatT angleAlpha = (3.0 * PI) / 4.0;

                     // interpolation coord
                     floatT interCoordx_1, interCoordy_1;
                     floatT interCoordx_2, interCoordy_2;
                     floatT interCoordy;

                     floatT distFromJetStart;

                     floatT Integral;

                     for(int step = 0; step < Steps; step++) {
                            // set integral to 0 for += notation
                            Integral = 0.0;

                            // calculate the slope from alpha
                     	floatT slope = std::tan(angleAlpha);

                            // point-slope form of equation for straight line y = m(x-x1)+y1
                            // calculate for each x value of the grid a value for y
                     	for(int interCoordx = _JetStartX; interCoordx >= 0; interCoordx--) {

                     		interCoordy = slope * (interCoordx - _JetStartX) + _JetStartY;

                                   // check wheather the point lies on the grid or not (needed for
                                   // integration which does not start at the origin (real (0.0,0.0))

                                   if(interCoordy < _smearEnerDensGrid.getMaxSitesPerDirection()
                                          && interCoordy >= 0.0) {
                                          // find the four gridpoints {(x1,y1),(x2,y1),(x1,y2),(x2,y2)}
                                          // for the interpolation around a point (x,y)
                                          interCoordx_1 = (int) std::floor((floatT) interCoordx);
                     			interCoordx_2 = (int) std::floor((floatT) interCoordx + 1.0);
                     			interCoordy_1 = (int) std::floor(interCoordy);
                     			interCoordy_2 = (int) std::floor(interCoordy + 1.0);

                                          // define the points where the energy density grid is read out
                     			Site siteP1(interCoordx_1, interCoordy_1);
                                          Site siteP2(interCoordx_1, interCoordy_2);
                                          Site siteP3(interCoordx_2, interCoordy_1);
                                          Site siteP4(interCoordx_2, interCoordy_2);
                                          Site siteP(interCoordx, interCoordy);

                     			floatT Interpol = _energyDensityInterpolation(
                                                 siteP1,
                                                 siteP2,
                                                 siteP3,
                                                 siteP4,
                                                 siteP);

                                          // calculate the distance from the jet origin
                     			distFromJetStart = std::hypot(
                                                 (((floatT) interCoordx - (floatT) _JetStartX) / 100.0),
                                                 ((interCoordy - (floatT) _JetStartY) / 100.0));
                                          // store point at line and the interpolated value at that point
                     			_RadiusTwo.at(_JetStartX - interCoordx) = distFromJetStart;
                     			_EDensTwo.at(_JetStartX - interCoordx) = Interpol;
                     		}
                     		else {
                                          // calculate the distance from the jet origin
                                          distFromJetStart = std::hypot(
                                                 (((floatT) interCoordx - (floatT) _JetStartX) / 100.0),
                                                 ((interCoordy - (floatT) _JetStartY) / 100.0));
                                          // store point at line and the interpolated value at that point
                                          _RadiusTwo.at(_JetStartX - interCoordx) = distFromJetStart;
                                          _EDensTwo.at(_JetStartX - interCoordx) = 0.0;
                     		}
                     	}

                            // perform a trapezoidal integration along the line
                            for(int inteStep = 0; inteStep < _JetStartX; inteStep++)
                            {
                                   Integral += 0.5*(_RadiusTwo.at(inteStep + 1) - _RadiusTwo.at(inteStep))
                                          * (_EDensTwo.at(inteStep) + _EDensTwo.at(inteStep + 1));
                            }

                     // calculate the angle from the center of mass (real (0.0,0.0))
                     interCoordy = slope * (0 -  _JetStartX) + (floatT) _JetStartY;

                     PhiTwo = std::abs(std::atan(
                            (((floatT) _smearEnerDensGrid.getMaxSitesPerDirection() - 1.)/2. - interCoordy) /
                            (((floatT) _smearEnerDensGrid.getMaxSitesPerDirection() - 1.)/2.)) + PI);

                     // store the values for the angle and the integral in an array
                     _AngleSec2.at(step) = PhiTwo;
                     _EDensSec2.at(step) = Integral;

                     angleAlpha += _AngleStep;
                     }
              }

              // calculation for the upper sector (sector3) from 1/4 pi to 3/4 pi
              void _sector3() {
                     // the roles of x and y must be changed because the slope/tan diverges at pi/2
                     // PhiThree is the angle with respect to the origin in real space
                     floatT PhiThree = 0.0;

                     // alpha is the angle seen from the start point of the integration
                     floatT angleAlpha = -PI / 4.0;

                     // interpolation coord
                     floatT interCoordx_1, interCoordy_1;
                     floatT interCoordx_2, interCoordy_2;
                     floatT interCoordx;

                     floatT distFromJetStart;

                     floatT Integral;

                     for(int step = 0; step < Steps; step++) {
                            // set integral to 0 for += notation
                            Integral = 0.0;

                            // calculate the slope from alpha
                     	floatT slope = std::tan(angleAlpha);

                            // point-slope form of equation for straight line x = m(y-y1)+x1
                            // calculate for each y value of the grid a value for x
                     	for(int interCoordy = _JetStartY;
                                    interCoordy < _smearEnerDensGrid.getMaxSitesPerDirection();
                                    interCoordy++) {

                                   floatT interCoordx = -slope * (interCoordy - _JetStartY)
                                                        + (floatT) _JetStartX;

                                   // check wheather the point lies on the grid or not (needed for
                                   // integration which does not start at the origin (real (0.0,0.0))
                     		if(interCoordx < _smearEnerDensGrid.getMaxSitesPerDirection()
                                          && interCoordx >= 0.0) {
                                          // find the four gridpoints {(x1,y1),(x2,y1),(x1,y2),(x2,y2)}
                                          // for the interpolation around a point (x,y)
                                          interCoordy_1 = (int) std::floor((floatT) interCoordy);
                     			interCoordy_2 = (int) std::floor((floatT) interCoordy + 1.0);
                     			interCoordx_1 = (int) std::floor(interCoordx);
                     			interCoordx_2 = (int) std::floor(interCoordx + 1.0);

                                          // define the points where the energy density grid is read out
                     			Site siteP1(interCoordx_1, interCoordy_1);
                                          Site siteP2(interCoordx_1, interCoordy_2);
                                          Site siteP3(interCoordx_2, interCoordy_1);
                                          Site siteP4(interCoordx_2, interCoordy_2);
                                          Site siteP(interCoordx, interCoordy);

                     			floatT Interpol = _energyDensityInterpolation(
                                                 siteP1,
                                                 siteP2,
                                                 siteP3,
                                                 siteP4,
                                                 siteP);

                                          // calculate the distance from the jet origin
                            		distFromJetStart = std::hypot(
                                                 ((interCoordx - (floatT) _JetStartX) / 100.0),
                                                 (((floatT)interCoordy - (floatT) _JetStartY) / 100.0));
                                          // store point at line and the interpolated value at that point
                     			_RadiusThree.at(interCoordy - _JetStartY) = distFromJetStart;
                     			_EDensThree.at(interCoordy - _JetStartY) = Interpol;
                     		}
                     		else {
                                          // calculate the distance from the jet origin
                            		distFromJetStart = std::hypot(
                                                 ((interCoordx - (floatT) _JetStartX) / 100.0),
                                                 (((floatT)interCoordy - (floatT) _JetStartY) / 100.0));
                                          // store point at line and the interpolated value at that point
                     			_RadiusThree.at(interCoordy - _JetStartY) = distFromJetStart;
                     			_EDensThree.at(interCoordy - _JetStartY) = 0.0;
                     		}
                     	}
                            // perform a trapezoidal integration along the line
                     	for(int inteStep = 0; inteStep < _smearEnerDensGrid.getMaxSitesPerDirection()
                                   - _JetStartY; inteStep++) {
                     		Integral += 0.5
                                          * (_RadiusThree.at(inteStep + 1) - _RadiusThree.at(inteStep))
                                          * (_EDensThree.at(inteStep) + _EDensThree.at(inteStep + 1));
                     	}

                     // calculate the angle from the center of mass (real (0.0,0.0))
                     interCoordx = slope * (_smearEnerDensGrid.getMaxSitesPerDirection()
                                   - _JetStartY) + (floatT) _JetStartX;
                     PhiThree = std::abs(std::atan(
                            (interCoordx - ((floatT) _smearEnerDensGrid.getMaxSitesPerDirection() - 1.)/2.) /
                             (((floatT) _smearEnerDensGrid.getMaxSitesPerDirection() - 1.)/2.)) + PI / 2.0);

                     // store the values for the angle and the integral in an array
                     _AngleSec3.at(step) = PhiThree;
                     _EDensSec3.at(step) = Integral;

                     angleAlpha += _AngleStep;
                     }
              }

              // calculation for the lower sector from 5/4 pi to 7/4 pi
              void _sector4() {
                     // the roles of x and y must be changed because the slope/tan diverges at pi/2
                     // PhiThree is the angle with respect to the origin in real space
                     floatT PhiFour = 0.0;

                     // alpha is the angle seen from the start point of the integration
                     floatT angleAlpha = -PI / 4.0;

                     // interpolation coord
                     floatT interCoordx_1, interCoordy_1;
                     floatT interCoordx_2, interCoordy_2;
                     floatT interCoordx;

                     floatT distFromJetStart;

                     floatT Integral;


                     for(int step = 0; step < Steps; step++) {
                            // set integral to 0 for += notation
                            Integral = 0.0;

                            // calculate the slope from alpha
                     	floatT slope = std::tan(angleAlpha);

                            //point-slope form of equation for straight line x = m(y-y1)+x1
                            //calculate for each y value of the grid a value for x
                     	for(int interCoordy = _JetStartY; interCoordy >= 0; interCoordy--)
                     	{
                     		floatT interCoordx = -slope * (interCoordy - _JetStartY)
                                                        + (floatT) _JetStartX;
                                   // check wheather the point lies on the grid or not (needed for
                                   // integration which does not start at the origin (real (0.0,0.0))
                     		if(interCoordx < _smearEnerDensGrid.getMaxSitesPerDirection()
                                          && interCoordx >= 0.0) {
                                          // find the four gridpoints {(x1,y1),(x2,y1),(x1,y2),(x2,y2)}
                                          // for the interpolation around a point (x,y)
                                          interCoordy_1 = (int) std::floor((floatT) interCoordy);
                     			interCoordy_2 = (int) std::floor((floatT) interCoordy + 1.0);
                     			interCoordx_1 = (int) std::floor(interCoordx);
                     			interCoordx_2 = (int) std::floor(interCoordx + 1.0);

                                          // define the points where the energy density grid is read out
                     			Site siteP1(interCoordx_1, interCoordy_1);
                                          Site siteP2(interCoordx_1, interCoordy_2);
                                          Site siteP3(interCoordx_2, interCoordy_1);
                                          Site siteP4(interCoordx_2, interCoordy_2);
                                          Site siteP(interCoordx, interCoordy);

                     			floatT Interpol = _energyDensityInterpolation(
                                                 siteP1,
                                                 siteP2,
                                                 siteP3,
                                                 siteP4,
                                                 siteP);

                                          // calculate the distance from the jet origin
                            		distFromJetStart = std::hypot(
                                                 ((interCoordx - (floatT) _JetStartX) / 100.0),
                                                 (((floatT) interCoordy - (floatT) _JetStartY) / 100.0));
                                          // store point at line and the interpolated value at that point

                     			_RadiusFour.at(_JetStartY - interCoordy) = distFromJetStart;
                     			_EDensFour.at(_JetStartY - interCoordy) = Interpol;
                     		}
                     		else {
                                          // calculate the distance from the jet origin
                                          distFromJetStart = std::hypot(
                                                 ((interCoordx - (floatT) _JetStartX) / 100.0),
                                                 (((floatT) interCoordy - (floatT) _JetStartY) / 100.0));
                                          // store point at line and the interpolated value at that point

                     			_RadiusFour.at(_JetStartY - interCoordy) = distFromJetStart;
                     			_EDensFour.at(_JetStartY - interCoordy) = 0.0;
                     		}
                     	}
                            //perform a trapezoidal integration along the line
                     	for(int inteStep = 0; inteStep < _JetStartY; inteStep++)
                     	{
                     		Integral += 0.5 * (_RadiusFour.at(inteStep + 1) - _RadiusFour.at(inteStep))
                                          * (_EDensFour.at(inteStep) + _EDensFour.at(inteStep + 1));
                     	}

                     //calculate the angle from the center of mass (real (0.0,0.0))
                     interCoordx = slope * (0 - _JetStartY) + (floatT) _JetStartX;
                     PhiFour = std::abs(std::atan(
                            (((floatT) _smearEnerDensGrid.getMaxSitesPerDirection() - 1.)/2. - interCoordx) /
                            (((floatT) _smearEnerDensGrid.getMaxSitesPerDirection() - 1.)/2.))
                            + (6.0 * PI) / 4.0);

                     // store the values for the angle and the integral in an array
                     _AngleSec4.at(step) = PhiFour;
                     _EDensSec4.at(step) = Integral;

                     angleAlpha += _AngleStep;
                     }
              }

              void _angles() {
                     for (int i = 0; i < Steps; i++) {
                            _AngleSec1.at(Steps+i) = _AngleSec3.at(i);
                     }
                     for (int i = 0; i < Steps; i++) {
                            _AngleSec1.at(2*Steps+i) = _AngleSec2.at(i);
                     }
                     for (int i = 0; i < Steps; i++) {
                            _AngleSec1.at(3*Steps+i) = _AngleSec4.at(i);
                     }
              }

              // bilinear interpolation function between the grid points in the energy density grid
              floatT inline _energyDensityInterpolation(
                                   Site intSiteLowerLeftCor,
                                   Site intSiteUpperLeftCor,
                                   Site intSiteLowerRightCor,
                                   Site intSiteUpperRightCor,
                                   Site intSite){

                     // compute distances between points of the square
                     floatT x2x = intSiteLowerRightCor.x() - intSite.x();
                     floatT x2x1 = intSiteLowerRightCor.x() - intSiteLowerLeftCor.x();
                     floatT xx1 = intSite.x() - intSiteLowerLeftCor.x();
                     floatT y2y = intSiteUpperRightCor.y() - intSite.y();
                     floatT y2y1 = intSiteUpperRightCor.y() - intSiteLowerRightCor.y();
                     floatT yy1 = intSite.y() - intSiteLowerRightCor.y();

                     // compute the interpolation value and return it
                     return  (_smearEnerDensGrid.getSite(intSiteLowerLeftCor) * x2x * y2y
                                   + _smearEnerDensGrid.getSite(intSiteLowerRightCor) * xx1 * y2y
                                   + _smearEnerDensGrid.getSite(intSiteUpperLeftCor) * x2x * yy1
                                   + _smearEnerDensGrid.getSite(intSiteUpperRightCor) * xx1 * yy1)
                                   /(x2x1 * y2y1);
              }

              floatT inline _integralNormalization() {
                     // normalization of the integral by summing the value of the integrals
                     // for all angles phi of one event
                     floatT IntegralNormalization = 0.0;

                     for(int i = 0; i < Steps; i++){
                            IntegralNormalization += ( _EDensSec1.at(i)
                                                     + _EDensSec2.at(i)
                                                     + _EDensSec3.at(i)
                                                     + _EDensSec4.at(i));

                     }
                     return IntegralNormalization;
              }

              void _averageIntegralsOverAllEvents(int NEvents) {
                     floatT IntegralNormalization = _integralNormalization();
                     for (int i = 0; i < Steps; i++) {
                            _Integral.at(i) += _EDensSec1.at(i) / (IntegralNormalization * NEvents);
                            _Integral.at(Steps + i) += _EDensSec3.at(i) / (IntegralNormalization * NEvents);
                            _Integral.at(2 * Steps + i) += _EDensSec2.at(i) / (IntegralNormalization * NEvents);
                            _Integral.at(3 * Steps + i) += _EDensSec4.at(i) / (IntegralNormalization * NEvents);
                     }
              }

       public:
              // constructor
              IntegratedEnergyDensity(Site startSite, EnergyDensity<floatT> & newEnergDens):
                                                        _RealX(startSite.x()),
                                                        _RealY(startSite.y()),
                                                        _JetStartX((int) ((startSite.x() + 15. ) * 100.) ),
                                                        _JetStartY((int) ((startSite.y() + 15. ) * 100.) ),
                                                        _AngleStep((PI / 2.0) / (Steps)),
                                                        _smearEnerDensGrid(
                                                               * newEnergDens.getSmearedEnergyDensData()),
                                                        _AngleSec1(4*Steps),
                                                        _EDensSec1(Steps),
                                                        _AngleSec2(Steps),
                                                        _EDensSec2(Steps),
                                                        _AngleSec3(Steps),
                                                        _EDensSec3(Steps),
                                                        _AngleSec4(Steps),
                                                        _EDensSec4(Steps),
                                                        _Integral(4*Steps),
                                                        _RadiusOne(newEnergDens.getSmearedEnergyDensData()
                                                                      -> getMaxSitesPerDirection()
                                                                      - startSite.x()),
                                                        _EDensOne(newEnergDens.getSmearedEnergyDensData()
                                                                      -> getMaxSitesPerDirection()
                                                                      - startSite.x()),
                                                        _RadiusTwo(std::abs(startSite.x()
                                                               - newEnergDens.getSmearedEnergyDensData()
                                                               -> getMaxSitesPerDirection())),
                                                        _EDensTwo(std::abs(startSite.x()
                                                               - newEnergDens.getSmearedEnergyDensData()
                                                               -> getMaxSitesPerDirection())),
                                                        _RadiusThree( newEnergDens.getSmearedEnergyDensData()
                                                               -> getMaxSitesPerDirection()
                                                               - startSite.y()),
                                                        _EDensThree( newEnergDens.getSmearedEnergyDensData()
                                                               -> getMaxSitesPerDirection()
                                                               - startSite.y()),
                                                        _RadiusFour(std::abs(startSite.y()
                                                               - newEnergDens.getSmearedEnergyDensData()
                                                               -> getMaxSitesPerDirection())),
                                                        _EDensFour(std::abs(startSite.y()
                                                               - newEnergDens.getSmearedEnergyDensData()
                                                               -> getMaxSitesPerDirection())) {}

              void computeIntegratedEnergyDensity(int NEvents){
              // compute the integrated energy density in each sector
              _sector1();
              _sector2();
              _sector3();
              _sector4();

              // normalize and merge the integral
              _averageIntegralsOverAllEvents(NEvents);

              // put computed angles of each sector in one output vector
              _angles();
       }

              std::vector<floatT> * getIntegrationEDensValAngles() {
              return & _AngleSec1;
       }

              std::vector<floatT> * getIntegrationEDensValIntegrals() {
                     return & _Integral;
              }

              friend class Eccentricity<floatT>;
              friend class FlowCoefficients<floatT>;
};

//TODO get rid of the sectors!!!
template<class floatT>
class Eccentricity{
       private:
              std::vector<std::complex<floatT>> _EccCounter;
              std::vector<std::complex<floatT>> _EccDenom;

              std::vector<floatT> _EventEccentricity;

              IntegratedEnergyDensity<floatT> & _intEnergDens;

              void _computeEccConterAndDnomPerSector(  std::vector<floatT> & radius,
                                          std::vector<floatT> & eneDens,
                                          std::vector<floatT> & angle ){
                     // define the imaginary I
                     const std::complex<floatT> I(0.0, 1.0);

                     // perform a trapezoidal integration to calculate the eccentricity e2 ... e5

                     for (int i = 0; i < Steps; i++) {
                            for(int n =  2; n <= 5; n++) {
                                   for(int j = 0; j < (int) radius.size()-1; j++) {
                                          _EccCounter.at(n - 2) += (floatT) 0.5 * (radius.at(j + 1) - radius.at(j))
                                                        * (floatT) (eneDens.at(j) * std::pow(radius.at(j),(n + 1))
                                                        + eneDens.at(j + 1)
                                                        * std::pow(radius.at(j + 1),(n + 1)))
                                                        * std::exp(I * ((floatT) n ) * angle.at(i));

                                          _EccDenom.at(n - 2) += (floatT) 0.5 * (radius.at(j + 1) - radius.at(j))
                                                        * (floatT) (eneDens.at(j) * std::pow(radius.at(j),(n + 1))
                                                        + eneDens.at(j + 1)
                                                        * std::pow(radius.at(j + 1),(n + 1)));
                                   }
                            }
                     }
              }

              void inline _computeEcc(){
                     // calculation of the eccentricity count and Dnorm for sector 1
                     _computeEccConterAndDnomPerSector(_intEnergDens._RadiusOne,  _intEnergDens._EDensOne,
                                    _intEnergDens._AngleSec1);

                     // calculation of the eccentricity count and Dnorm for sector 2
                     _computeEccConterAndDnomPerSector(_intEnergDens._RadiusTwo, _intEnergDens._EDensTwo,
                                   _intEnergDens._AngleSec2);

                     // calculation of the eccentricity count and Dnorm for sector 3
                     _computeEccConterAndDnomPerSector(_intEnergDens._RadiusThree, _intEnergDens._EDensThree,
                                   _intEnergDens._AngleSec3);

                     // calculation of the eccentricity count and Dnorm for sector 4
                     _computeEccConterAndDnomPerSector(_intEnergDens._RadiusFour, _intEnergDens._EDensFour,
                                   _intEnergDens._AngleSec4);

                     // calculation of the eccentricity
                     for(int n = 0; n < 4; n++){
                            _EventEccentricity.at(n) = std::abs(-_EccCounter.at(n)/_EccDenom.at(n));
                            _EccCounter.at(n) = 0;
                            _EccDenom.at(n) = 0;
                     }
              }

       public:
              // constructor
              Eccentricity(IntegratedEnergyDensity<floatT> & newIntEnergDens):
                            _EccCounter(4),
                            _EccDenom(4),
                            _EventEccentricity(4),
                            _intEnergDens(newIntEnergDens){}

              // calculate the eccentricity e2 ... e5 for one event and store the value in the array
              void computeEccentricity(){
                     _computeEcc();
              }

              std::vector<floatT> * getEccenetricityData(){
                     return & _EventEccentricity;
              }
};

template<class floatT>
class FlowCoefficients{
       private:
              std::vector<floatT> _MergedAngles;
              std::vector<floatT> _MergedIntegrals;

              std::vector<floatT> _Fourier;

              IntegratedEnergyDensity<floatT> & _intEnergDens;

              // merged vectors
              // the order of the sectors is adapted to a mathematical positive integration direction
              // sector 1 --> sector3 --> sector2 --> sector 4
              void _mergeVectors() {
                     for(int i = 0; i < Steps; i++) {
                            _MergedAngles.at(i) = _intEnergDens._AngleSec1.at(i);
                     }
                     for(int i = 0; i < Steps; i++) {
                            _MergedAngles.at(Steps + i) = _intEnergDens._AngleSec3.at(i);
                     }
                     for(int i = 0; i < Steps; i++) {
                            _MergedAngles.at(2 * Steps + i) = _intEnergDens._AngleSec2.at(i);
                     }
                     for(int i = 0; i < Steps; i++) {
                            _MergedAngles.at(3 * Steps + i) = _intEnergDens._AngleSec4.at(i);
                     }

                     for(int i = 0; i < Steps; i++) {
                            _MergedIntegrals.at(i) = _intEnergDens._EDensSec1.at(i);
                     }
                     for(int i = 0; i < Steps; i++) {
                            _MergedIntegrals.at(Steps + i) = _intEnergDens._EDensSec3.at(i);
                     }
                     for(int i = 0; i < Steps; i++) {
                            _MergedIntegrals.at(2 * Steps + i) = _intEnergDens._EDensSec2.at(i);
                     }
                     for(int i = 0; i < Steps; i++) {
                            _MergedIntegrals.at(3 * Steps + i) = _intEnergDens._EDensSec4.at(i);
                     }
              }

              //calculation of the flow coefficients
              void _flowCoefficients() {
                     // Compute an Fourier series
                     floatT fourCoeff_0 = 0.0;
                     for(int i=0; i < 4 * Steps - 1 ; i++) {
                            fourCoeff_0 += (
                                   ((_MergedIntegrals.at(i) + _MergedIntegrals.at(i + 1)) / 2.0)
                                   * _intEnergDens._AngleStep) / PI;
                     }

                     floatT aFourier = 0.0;
                     for(int aCount=1; aCount <= 3; aCount++) {
                            aFourier = 0.0;
                            for(int i=0; i < 4 * Steps - 1; i++) {
                                   aFourier += (
                                          (((_MergedIntegrals.at(i) + _MergedIntegrals.at(i + 1)) / 2.0)
                                                 * std::cos((floatT) aCount
                                                 * ((_MergedAngles.at(i) + _MergedAngles.at(i + 1)) / 2.0)))
                                          * _intEnergDens._AngleStep) / PI;
                            }
                            _Fourier.at(aCount - 1) = aFourier / (fourCoeff_0 / 2.0);
                            aFourier = 0.0;
                     }

                     floatT bFourier = 0.0;
                     for(int bCount=1; bCount <= 3; bCount++) {
                            bFourier = 0.0;
                            for(int k=0; k < 4*Steps-1; k++) {
                                   bFourier += ((((_MergedIntegrals.at(k)
                                                        + _MergedIntegrals.at(k + 1)) / 2.0)
                                                        * std::sin((floatT) bCount
                                                        * ((_MergedAngles.at(k)
                                                        + _MergedAngles.at(k + 1)) / 2.0)))
                                                        * _intEnergDens._AngleStep) / PI;
                            }
                            _Fourier.at(bCount + 2) = bFourier / (fourCoeff_0 / 2.0);
                            bFourier = 0.0;
                     }
                     fourCoeff_0 = 0.0;
              }

       public:
              // constructor
              FlowCoefficients(IntegratedEnergyDensity<floatT> & newIntEnergDens):
                            _MergedAngles(4 * Steps),
                            _MergedIntegrals(4 * Steps),
                            _Fourier(6),
                            _intEnergDens(newIntEnergDens) {}

              void computeFlowCoefficients(){
                     // merge them
                     _mergeVectors();

                     // compute flow coefficients
                     _flowCoefficients();
              }

              std::vector<floatT> * getFlowCoefficients(){
                     return & _Fourier;
              }
};

void computeMultipleEvents( Grid<PREC> & rawDataGrid,
                            FileReader<PREC> & fr,
                            EnergyDensity<PREC> & energDens,
                            IntegratedEnergyDensity<PREC> & intEnergDens,
                            Eccentricity<PREC> & ecc,
                            FlowCoefficients<PREC> & flCoeffs,
                            int Event,
                            int NEvents,
                            int NNucleonsCore,
                            PREC NNCross) {

       energDens.energyDensity(rawDataGrid);

       energDens.smearedEnergyDensity();

       intEnergDens.computeIntegratedEnergyDensity(NEvents);

       ecc.computeEccentricity();

       flCoeffs.computeFlowCoefficients();

}

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

       // define size of grid = elems*elems + 2 elems
       int elems = 3001;

       // define a start site for the integrated energy density computation
       Site startSite(0,0);

       FileReader<PREC> fr(NEvents, NNucleonsCore, "Pb67.6.txt");
       fr.readData(NEvents, NNucleonsCore);

       std::string filename3 = "Eccentricities.dat";
       FileWriter<PREC> file3(filename3);

       std::string filename4 = "FourierCoeff.dat";
       FileWriter<PREC> file4(filename4);

       std::vector<PREC> globalIntegral(4*Steps);
       std::vector<PREC> globalAngle(4*Steps);

       #pragma omp parallel for
       for (int Event = 1; Event <= NEvents; Event++) {
              std::cout << Event << '\n';
              Grid<PREC> rawDataGrid(elems);
              EnergyDensity<PREC> energDens(NNCross, elems);
              IntegratedEnergyDensity<PREC> intEnergDens(startSite, energDens);
              Eccentricity<PREC> ecc(intEnergDens);
              FlowCoefficients<PREC> flCoeffs(intEnergDens);
              #pragma omp critical
              fr.getData(Event, NNucleonsCore, rawDataGrid);
              computeMultipleEvents(rawDataGrid, fr, energDens, intEnergDens, ecc, flCoeffs, Event, NEvents, NNucleonsCore, NNCross);
              #pragma omp critical
              for (int i = 0; i < 4*Steps; i++) {
                     globalIntegral.at(i) += intEnergDens.getIntegrationEDensValIntegrals() -> at(i);
              }
              // write eccentricities in "Eccentricities.dat"
              file3 << Event << "\t"  << ecc.getEccenetricityData() -> at(0) << "\t"
                     << ecc.getEccenetricityData() -> at(1) << "\t"
                     << ecc.getEccenetricityData() -> at(2) << "\t"
                     << ecc.getEccenetricityData() -> at(3) << std::endl;
              // write the flow coefficients in "FourierCoef.dat";
              file4  << flCoeffs.getFlowCoefficients() -> at(0) << "\t"
                     << flCoeffs.getFlowCoefficients() -> at(1) << "\t"
                     << flCoeffs.getFlowCoefficients() -> at(2) << "\t"
                     << flCoeffs.getFlowCoefficients() -> at(3) << "\t"
                     << flCoeffs.getFlowCoefficients() -> at(4) << "\t"
                     << flCoeffs.getFlowCoefficients() -> at(5) << std::endl;
              // write one example smeared energy density
              if(Event == 1){
                     std::string filename1 = "EDensity.dat";
                     FileWriter<PREC> file1(filename1);
                     // write smeared Energy density in "EDensity.dat"
                     Site site(0,0);
                     do{
                            file1  << site.x()/100. - 15. << "\t"
                                   << site.y()/100. - 15. << "\t"
                                   << energDens.getSmearedEnergyDensData() -> getSite(site) << std::endl;
                     }while(energDens.getSmearedEnergyDensData() -> runThroughGrid(site));
                     globalAngle = * intEnergDens.getIntegrationEDensValAngles();
              }
       }

       std::string filename2 = "Integration.dat";
       FileWriter<PREC> file2(filename2);

       // write angles and integrated energy density in "Integration.dat"
       for (int i = 0; i < 4 * Steps; i++) {
              file2  << globalAngle.at(i) << "\t"
                     << globalIntegral.at(i) << std::endl;
       }

       return 0;
}



// Testing routins
// Test indexer:
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

//Test runThroughGrid()
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

// Test FileWriter
void testFileWriter(int elems, PREC NNCross, int NNucleonsCore, int NEvents){
       Site startSite(0,0);

       FileReader<PREC> fr(NEvents, NNucleonsCore, "Pb67.6.txt");
       fr.readData(NEvents, NNucleonsCore);

       std::string test = "T1.dat";
       FileWriter<PREC> test1(test);

              #pragma omp parallel for
              for (int Event = 1; Event <= NEvents; Event++) {
                     std::cout << Event << '\n';
                     Grid<PREC> rawDataGrid(elems);
                     fr.getData(Event, NNucleonsCore, rawDataGrid);

                     EnergyDensity<PREC> energDens(NNCross, elems);
                     IntegratedEnergyDensity<PREC> intEnergDens(startSite, energDens);
                     Eccentricity<PREC> ecc(intEnergDens);

                     energDens.energyDensity(rawDataGrid);

                     energDens.smearedEnergyDensity();

                     intEnergDens.computeIntegratedEnergyDensity(NEvents);

                     ecc.computeEccentricity();

                     test1 << Event << "\t" << ecc.getEccenetricityData() -> at(0) << "\t"
                            << ecc.getEccenetricityData() -> at(1) << "\t"
                            << ecc.getEccenetricityData() -> at(2) << "\t"
                            << ecc.getEccenetricityData() -> at(3) << std::endl;
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

// Test Integrated Energy Density
void testIntegratedEnergyDensity(){
       int elems = 3001;

       //define some parameter which are necessary for the computation. For the testing purpose they can be random.
       PREC NNCross = 67.6;
       int NEvents = 102;
       int NNucleonsCore = 208;

       Site site(0,0);

       //you can varie this for testing purpose
       Site startSite(0,0);

       //fill the grid with values
       FileReader<PREC> fr(NEvents, NNucleonsCore, "Pb67.6.txt");
       fr.readData(NEvents, NNucleonsCore);

       //do the test often
       //#pragma omp parallel for
       for(int event = 100; event <= NEvents; event++){
              Grid<PREC> grid(elems);
              std::cout << event << '\n';

              fr.getData(event, NNucleonsCore, grid);

              Site si(0,0);
              do{
                     if(grid.getSite(si)  != 0){
                            si.printSite();
                            std::cout << grid.getSite(si) << '\n';
                     }
              }while(grid.runThroughGrid(si));

              //the Integrated energy density needs a energy density which needs to have no bugs, tested before!
              EnergyDensity<PREC> energDens(NNCross, elems);

              energDens.energyDensity(grid);

              energDens.smearedEnergyDensity();

              IntegratedEnergyDensity<PREC> intEnergDens(startSite, energDens);

              intEnergDens.computeIntegratedEnergyDensity(NEvents);

              Eccentricity<PREC> ecc(intEnergDens);

              ecc.computeEccentricity();

              for (size_t i = 0; i < 4*Steps; i++) {
                     if((*intEnergDens.getIntegrationEDensValIntegrals())[i]
                            != (*intEnergDens.getIntegrationEDensValIntegrals())[i]){
                     //std::cout << i << "\t" << (*intEnergDens.getIntegrationEDensValIntegrals())[i] << '\n';
                     }
              }

              if(ecc.getEccenetricityData() -> at(0) != ecc.getEccenetricityData() -> at(0)){
              //       std::cout << ecc.getEccenetricityData() -> at(0);
              }
              if(ecc.getEccenetricityData() -> at(1) != ecc.getEccenetricityData() -> at(1)){
              //       std::cout << ecc.getEccenetricityData() -> at(1);
              }
              if(ecc.getEccenetricityData() -> at(2) != ecc.getEccenetricityData() -> at(2)){
              //       std::cout << ecc.getEccenetricityData() -> at(2);
              }
              if(ecc.getEccenetricityData() -> at(3) != ecc.getEccenetricityData() -> at(3)){
              //       std::cout << ecc.getEccenetricityData() -> at(3) << '\n';
              }
       }
}
