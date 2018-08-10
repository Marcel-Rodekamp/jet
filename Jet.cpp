#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<cmath>
#include<complex>
#include<vector>
#include<ctime>

#define PREC double
#define PI 3.14159265358979323846
#define Steps 500

class Site;
template<class floatT> class Grid;
template<class floatT> class FileWriter;
template<class floatT> class EnergyDensity;
template<class floatT> class IntegratedEnergyDensity;
template<class floatT> class Eccentricity;
template<class floatT> class FlowCoefficients;


void indexerTest(int max);
void runThroughGridTest(int max);
void computeTest();
void FileWriterTest();

// structure to define the lattice points leaving out the z component due to the contemplated event plane
class Site{
       private:
              std::vector<int> VCoordinates;

       public:
              // constructor
              Site(int x, int y){
                     VCoordinates.push_back(x);
                     VCoordinates.push_back(y);
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
              // TODO What happens at site s(0,0)!!!
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
                            if(newCoordY < _maxSitesPerDirection){
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
                                 std::cout << "File could not be opened" << '\n';
                                 return;
                          }

                         // set high precision
                         _fileStream.precision(15);
                         _fileStream.setf(std::ios::scientific);
                  }

       public:

              // constructor
              FileWriter(std::string fname) :
                            _fileName(fname) {
                     init();
              };

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

                     while(_grid.runThroughGrid(site)){
                            NColl = (int) grid.getSite(site);

                            // (15*800/46^4) * NColl^4
                            EDens = 0.002680093338717343 * NColl * NColl * NColl * NColl;

                            _grid.setSite(site,EDens);
                     }
              }

              void smearedEnergyDensity(){
                     Site site(0, 0);
                     Site tmpSite(0,0);
                     const int tmpRadNucleon = (int) _RadNucleon;

                     // test site (0,0)
                     if (_grid.getSite(site) != 0 ) {

                            // loop through a square with side length tmpRadNucleon centered on the
                            // not zero site
                            for(int x = site.x() - tmpRadNucleon; x < site.x() + tmpRadNucleon; x++){
                            for(int y = site.y() - tmpRadNucleon; y < site.y() + tmpRadNucleon; y++){
                                   tmpSite.setX(x);
                                   tmpSite.setY(y);
                                   _gridSmeared.addSite(tmpSite,_grid.getSite(site));
                            }
                            }
                     }

                     while(_grid.runThroughGrid(site)){
                            if (_grid.getSite(site) != 0 ) {

                                   // loop through a square with side length tmpRadNucleon centered on the
                                   // not zero site
                                   for(int x = site.x() - tmpRadNucleon; x < site.x() + tmpRadNucleon; x++){
                                   for(int y = site.y() - tmpRadNucleon; y < site.y() + tmpRadNucleon; y++){
                                          tmpSite.setX(x);
                                          tmpSite.setY(y);
                                          _gridSmeared.addSite(tmpSite,_grid.getSite(site));
                                   }
                                   }
                            }
                     }
              }

              Grid<floatT> * getSmearedEnergyDensData(){
                     return & _gridSmeared;
              }

              friend class IntegratedEnergyDensity<floatT>;

};

//TODO _realX usw in eine Site verwandeln
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

              // vectors for the averaged energy density
              std::vector<floatT> _AverageEDens1;
              std::vector<floatT> _AverageEDens2;
              std::vector<floatT> _AverageEDens3;
              std::vector<floatT> _AverageEDens4;

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
                                    interCoordx < _smearEnerDensGrid.getMaxSitesPerDirection();
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

                     PhiOne = std::abs(std::atan(
                            (interCoordy - ((floatT) _smearEnerDensGrid.getMaxSitesPerDirection() - 1.)/2.)
                                   /(((floatT) _smearEnerDensGrid.getMaxSitesPerDirection() - 1.)/2.)));

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

              void _normalizeIntegral() {
                     floatT IntegralNormalization = _integralNormalization();
                     // divide each integral by the normalization and the number of events
                     for(int j = 0; j < Steps; j++){
                            _AverageEDens1.at(j) += _EDensSec1.at(j) / IntegralNormalization;
                            _AverageEDens2.at(j) += _EDensSec2.at(j) / IntegralNormalization;
                            _AverageEDens3.at(j) += _EDensSec3.at(j) / IntegralNormalization;
                            _AverageEDens4.at(j) += _EDensSec4.at(j) / IntegralNormalization;
                     }
              }

              void _angles() {
                     for (int i = 0; i < Steps; i++) {
                            _AngleSec1.push_back(_AngleSec2.at(i));
                            _AngleSec1.push_back(_AngleSec3.at(i));
                            _AngleSec1.push_back(_AngleSec4.at(i));
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
                     return  (intSiteLowerLeftCor * x2x * y2y
                                   + intSiteLowerRightCor * xx1 * y2y
                                   + intSiteUpperLeftCor * x2x * yy1
                                   + intSiteUpperRightCor * xx1 * yy1)
                                   /(x2x1 * y2y1);
              }

              floatT inline _integralNormalization() {
                     // normalization of the integral by summing the value of the integrals
                     // for all angles phi of one event
                     floatT IntegralNormalization = 0.;

                     for(int i = 0; i < Steps; i++){
                            IntegralNormalization += (  _EDensSec1.at(i)
                                                     + _EDensSec2.at(i)
                                                     + _EDensSec3.at(i)
                                                     + _EDensSec4.at(i));
                     }

                     return IntegralNormalization;
              }

       public:
              // constructor
              IntegratedEnergyDensity(Site startSite, EnergyDensity<floatT> & newEnergDens):
                                                        _RealX(startSite.x()),
                                                        _RealY(startSite.y()),
                                                        _JetStartX((int) ((startSite.x() + 15. ) * 100.) ),
                                                        _JetStartY((int) ((startSite.y() + 15. ) * 100.) ),
                                                        _AngleStep((PI / 2.0) / (Steps)),
                                                        _smearEnerDensGrid(newEnergDens),
                                                        _AngleSec1(Steps),
                                                        _EDensSec1(Steps),
                                                        _AngleSec2(Steps),
                                                        _EDensSec2(Steps),
                                                        _AngleSec3(Steps),
                                                        _EDensSec3(Steps),
                                                        _AngleSec3(Steps),
                                                        _EDensSec3(Steps),
                                                        _Integral(4*Steps),
                                                        _AverageEDens1(Steps),
                                                        _AverageEDens2(Steps),
                                                        _AverageEDens3(Steps),
                                                        _AverageEDens4(Steps),
                                                        _RadiusOne(3000 - (int) ((startSite.x() + 15. ) * 100.)),
                                                        _EDensOne(3000 - (int) ((startSite.x() + 15. ) * 100.)),
                                                        _RadiusTwo((int) ((startSite.x() + 15. ) * 100.)),
                                                        _EDensTwo((int) ((startSite.x() + 15. ) * 100.)),
                                                        _RadiusThree(3000 - (int) ((startSite.y() + 15. ) * 100.)),
                                                        _EDensThree(3000 - (int) ((startSite.y() + 15. ) * 100.)),
                                                        _RadiusFour((int) ((startSite.y() + 15. ) * 100.)),
                                                        _EDensFour((int) ((startSite.y() + 15. ) * 100.)) {}

       void computeIntegratedEnergyDensity(){
              // compute the integrated energy density in each sector
              _sector1();
              _sector2();
              _sector3();
              _sector4();

              // normalize and merge the integral
              _normalizeIntegral();

              // put computed angles of each sector in one output vector
              _angles();
       }

       void averageIntegralsOverAllEvents(int NEvents) {
              for (int i = 0; i < Steps; i++) {
                     _Integral.at(i) = _AverageEDens1.at(i) / NEvents;
                     _Integral.at(Steps + i) = _AverageEDens2.at(i) / NEvents;
                     _Integral.at(2 * Steps + i) = _AverageEDens3.at(i) / NEvents;
                     _Integral.at(3 * Steps + i) = _AverageEDens4.at(i) / NEvents;
              }
       }

       std::vector<floatT> * getIntegrationEDensValAngles() {
              return & _AngleSec1;
       }

       std::vector<floatT> * getIntegrationEDensValIntegrals() {
              return & _Integral;
       }


              void angles() {
                     for (int i = 0; i < Steps; i++) {
                            _Angle.at(i) = _AngleSec1.at(i);
                            _Angle.at(Steps + i) = _AngleSec2.at(i);
                            _Angle.at(2 * Steps + i) = _AngleSec3.at(i);
                            _Angle.at(3 * Steps + i) = _AngleSec4.at(i);
                     }
              }

              void integrals(int NEvents) {
                     for (int i = 0; i < Steps; i++) {
                            _Integral.at(i) = _AverageEDens1.at(i) / NEvents;
                            _Integral.at(Steps + i) = _AverageEDens2.at(i) / NEvents;
                            _Integral.at(2 * Steps + i) = _AverageEDens3.at(i) / NEvents;
                            _Integral.at(3 * Steps + i) = _AverageEDens4.at(i) / NEvents;
                     }
              }

              std::vector<floatT> * getIntegrationEDensValAngles() {
                     return & _Angle;
              }

              std::vector<floatT> * getIntegrationEDensValIntegrals() {
                     return & _Integral;
              }

              friend class Eccentricity<floatT>;
              friend class FlowCoefficients<floatT>;
};

template<class floatT>
class Eccentricity{
       private:
              std::vector<std::complex<floatT>> _EccCounter;
              std::vector<std::complex<floatT>> _EccDenom;

              std::vector<floatT> _EventEccentricity;

              IntegratedEnergyDensity<floatT> & _intEnergDens;

              void _computeEccConterAndDnomPerSector(  std::vector<floatT> * radius,
                                          std::vector<floatT> * eneDens,
                                          std::vector<floatT> * angle ){
                     // define the imaginary I
                     const std::complex<floatT> I(0.0,1.0);

                     // perform a trapezoidal integration to calculate the eccentricity e2 ... e5
                     for (int i = 0; i < Steps; i++) {
                            for(int n =  2; n <= 5; n++) {
                                   for(int j = 0; j < 3000 - _intEnergDens._JetStartX; j++) {
                                          _EccCounter.at(n - 2) += 0.5 * (radius.at(j + 1) - radius.at(j))
                                                        * (eneDens.at(j) * std::pow(radius.at(j),(n + 1))
                                                               + eneDens.at(j + 1)
                                                               * std::pow(radius.at(j + 1),(n + 1)))
                                                        * exp(I * ((floatT) n * angle.at(i)));

                                          _EccDenom.at(n - 2) += 0.5 * (radius.at(j + 1) - radius.at(j))
                                                        * (eneDens.at(j)
                                                               * std::pow(radius.at(j),(n + 1))
                                                               + eneDens.at(j + 1)
                                                               * std::pow(radius.at(j + 1),(n + 1)));
                                   }
                            }
                     }
              }

              void inline _computeEcc(){
                     // calculation of the eccentricity count and Dnorm for sector 1
                     _computeEccConterAndDnomPerSector(_intEnergDens._RadiusOne, _intEnergDens._EDensOne,
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
                     }
              }

       public:
              // constructor
              Eccentricity(IntegratedEnergyDensity<floatT> & newIntEnergDens):
                            _intEnergDens(newIntEnergDens),
                            _EccCounter(4),
                            _EccDenom(4),
                            _EventEccentricity(4) {}

              // calculate the eccentricity e2 ... e5 for one event and store the value in the array
              void computeEccentricity(){
                     _computeEcc();
              }

              void eccentricitySector4() {
                     //define the imaginary I.
                     const std::complex<floatT> I(0.0,1.0);
                     // perform a trapezoidal integration to calculate the eccentricity e2 ... e5
                     for (int i = 0; i < Steps; i++) {
                            for(int n =  2; n <= 5; n++) {
                                   for(int j = 0; j < 3000 - _intEnergDens._JetStartX; j++) {
                                          _EccCounter4.at(n - 2) += (_intEnergDens._RadiusFour.at(j + 1)
                                                 - _intEnergDens._RadiusFour.at(j)) * 0.5
                                                 * (_intEnergDens._EDensFour.at(j)
                                                 * std::pow(_intEnergDens._RadiusFour.at(j),(n + 1))
                                                 + _intEnergDens._EDensFour.at(j + 1)
                                                 * std::pow(_intEnergDens._RadiusFour.at(j + 1),(n + 1)))
                                                 * exp(I * ((floatT) n * _intEnergDens._AngleSec4.at(i)));

                                          _EccDenom4.at(n - 2) += (_intEnergDens._RadiusFour.at(j + 1)
                                                 - _intEnergDens._RadiusFour.at(j)) * 0.5
                                                 * (_intEnergDens._EDensFour.at(j)
                                                 * std::pow(_intEnergDens._RadiusFour.at(j),(n + 1))
                                                 + _intEnergDens._EDensFour.at(j + 1)
                                                 * std::pow(_intEnergDens._RadiusFour.at(j + 1),(n + 1)));
                                   }
                            }
                     }
              }

              void eventEccentricity() {
                     //calculate the eccentricity e2 ... e5 for one event and store the value in the array
                     for(int n = 0; n < 4; n++){
                                  _EventEccentricity.at(n) = std::abs(
                                          -(_EccCounter1.at(n) + _EccCounter2.at(n)
                                          + _EccCounter3.at(n) + _EccCounter4.at(n))
                                          / (_EccDenom1.at(n) + _EccDenom2.at(n)
                                          + _EccDenom3.at(n) + _EccDenom4.at(n)));
                     }
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
              void _mergeVectors() {
                     for(int i = 0; i < Steps; i++) {
                            _MergedAngles.at(i) = _intEnergDens._AngleSec1.at(i);
                            _MergedAngles.at(Steps + i) = _intEnergDens._AngleSec2.at(i);
                            _MergedAngles.at(2 * Steps + i) = _intEnergDens._AngleSec3.at(i);
                            _MergedAngles.at(3 * Steps + i) = _intEnergDens._AngleSec4.at(i);

                            _MergedIntegrals.at(i) = _intEnergDens._EDensSec1.at(i);
                            _MergedIntegrals.at(Steps + i) = _intEnergDens._EDensSec2.at(i);
                            _MergedIntegrals.at(2 * Steps + i) = _intEnergDens._EDensSec3.at(i);
                            _MergedIntegrals.at(3 * Steps + i) = _intEnergDens._EDensSec4.at(i);
                     }
              }

              //calculation of the flow coefficients
              void _flowCoefficients() {
                     // Compute an Fourier series
                     floatT fourCoeff_0 = 0.0;
                     for(int i=0; i < 4 * Steps - 1; i++) {
                            fourCoeff_0 += (
                                   ((_MergedIntegrals.at(i)
                                          + _MergedIntegrals.at(i + 1)) / 2.0)
                                   * _intEnergDens._AngleStep) / PI;
                     }

                     floatT aFourier = 0.0;
                     for(int aCount=1; aCount <= 3; aCount++) {
                            aFourier = 0.0;
                            for(int i=0; i < 4 * Steps - 1; i++) {
                                   aFourier += (
                                          (((_MergedIntegrals.at(i) + _MergedIntegrals.at(i + 1)) / 2.0)
                                                 * std::cos((floatT) aCount
                                                 * ((_MergedAngles.at(i)+ _MergedAngles.at(i + 1)) / 2.0)))
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
              FlowCoefficients(EnergyDensity<floatT> & newEnergDens):
                            _intEnergDens(newEnergDens),
                            _MergedAngles(4 * Steps),
                            _MergedIntegrals(4 * Steps),
                            _Fourier(6) {}

              void computFlowCoefficients(){
                     // compute flow coefficients
                     _flowCoefficients();

                     // merge them
                     _mergeVectors();
              }

};


template<class floatT>
void readData(Grid<floatT> & grid, const std::string fileName, int NEvents, int NNucleonsCore ){
       std::fstream file;
       file.open(fileName.c_str(), std::ios::in);

       if(!file.is_open()){return;}

       std::vector<floatT> row(4);
       // loop through the data file which format is clarified by
       // x corrd \t y coord \t z coord \t NColl
       for(int i = 0; i < 2 * NEvents * NNucleonsCore; ++i){
              for(int col = 0; col < 4 ; ++col){
                     file >> row[col];
              }
              // translate the euclidean coordinates into the grid coordinates
              // leaving out the z component
              int x = (row[0] + 15)*100;
              int y = (row[1] + 15)*100;
              if(x > grid.getMaxSitesPerDirection() || y > grid.getMaxSitesPerDirection()) {
                     std::cout << "ERROR@readFile: Coordinates out of grid!" << '\n';
              }
              // compute the site
              Site site(x, y);
              // set NColl on the grid
              grid.setSite(site, row[3]);
       }

        file.close();
}

int main(int argc, char const *argv[]) {
       std::clock_t start;

       start = std::clock();


       // runThroughGridTest(3);

       // computeTest();

       // FileWriterTest();

       std::cout << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << "s" << '\n';

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

//Test runThroughGrid()

void runThroughGridTest(int max){
       Grid<PREC> grid(max);

       std::cout << "runThroughGrid(Site)" << '\n';

       Site site(0,0);
       int vectorIndex = 0;

       while(grid.runThroughGrid(site)){
              std::cout << "now at (x,y): (" << site.x() << "," << site.y() << ")" << '\n';
       }

       std::cout << "\nloop through lattice\n" << '\n';

       for (size_t y = 0; y <= max; y++) {
              for (size_t x = 0; x <= max; x++){
                     Site si(x,y);
                     std::cout << "now at (x,y): (" << si.x() << "," << si.y() << ")" << '\n';
              }
       }
}

// Test Energy Density computation and runtime

void computeTest(){
       // number of events
       int NEvents = 1;

       // std::cout << "Number of events: \n";

       // std::cin >> NEvents;

       // number of nucleons in the core of the element, e.g. 208 for Pb
       int NNucleonsCore = 208;

       // std::cout << "Number of nucleons per core: \n";

       // std::cin >> NNucleonsCore;

       // nucleon nucleon cross section for nucleon radius
       PREC NNCross = 67.6;

       // std::cout << "Nucleon Nucleon cross section in mb: \n";

       // std::cin >> NNCross;

       int elems = 3001;

       std::cout << "Initializing the data grid" << '\n';

       Grid<PREC> rawDataGrid(elems);

       std::cout << "Initializing the energy density computation" << '\n';

       EnergyDensity<PREC> energDens(NNCross, rawDataGrid.getMaxSitesPerDirection());

       std::cout << "Read data " << '\n';

       std::string filename = "Pb67.6.txt";

       FileWriter<PREC> file(filename);

       // file.readFile(rawDataGrid);

       std::cout << "Compute energy density" << '\n';

       energDens.energyDensity(rawDataGrid);

       std::cout << "Smear energy density" << '\n';

       energDens.smearedEnergyDensity();

       file.writeFileGrid(energDens.getSmearedEnergyDensData());

       std::cout << "Compute integrations in all directions" << '\n';

       Site startSite(0,0);

       IntegratedEnergyDensity<PREC> intEnergDens(startSite, energDens);



       file.writeFileVector(intEnergDens.getIntegrationEDensValAngles(), "Angle.dat");
       file.writeFileVector(intEnergDens.getIntegrationEDensValIntegrals(), "IntegratedEnergyDensity.dat");

       /*       Eccentricity<PREC> ecc(intEnergDens);

       ecc.eccentricitySector1();
       ecc.eccentricitySector2();
       ecc.eccentricitySector3();
       ecc.eccentricitySector4();
       ecc.eventEccentricity();

       file.writeFileVector(, "Eccentricity.dat");*/

       // file.writeFileGrid(energDens.getSmearedEnergyDensData());

       // file.writeFileGrid(energDens.getSmearedEnergyDensData());
}

// Test FileWriter

void FileWriterTest(){

       FileWriter<PREC> file("test.txt");

       file << "Dies is ein Test" << std::endl;

       Grid<PREC> grid(3001);

       readData<PREC>(grid, "Pb67.6.txt", 1, 208);
}
