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

struct Site;
template<class floatT> class Grid;

void indexerTest(int max);
void runThroughGridTest(int max);
void computeTest();

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

              void writeFileGrid(Grid<floatT> * grid, int accuracy = 10){
                     std::fstream file;
                     std::string newFilename;
                     newFilename.append("smearedEnergyDensity_");
                     newFilename.append(_filename);
                     file.open(newFilename, std::ios::out);

                     if(!file.is_open()){return;}

                     for(int x = 0; x < grid -> getMaxSitesPerDirection(); x += accuracy){
                            for(int y = 0; y < grid -> getMaxSitesPerDirection(); y += accuracy){

                            Site site(x,y);
                            file << ((floatT) x / 100.0) - 15.0 << "\t" << ((floatT) y / 100.0) - 15.0
                            << "\t" << grid -> getSite(site) << std::endl;


                            }
                     }

                     file.close();
              }

              void writeFileVector( std::vector<floatT> * data, std::string additionalFilename){
                     //open the file for jet data
                     std::fstream file;
                     std::string newFilename;
                     newFilename.append(additionalFilename);
                     newFilename.append(_filename);

                     for(int dataVectorIndex = 0; dataVectorIndex < data -> size(); dataVectorIndex++){
                            file << data -> at(dataVectorIndex) << std::endl;
                     }

                     //close the jet data file
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

};

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

              EnergyDensity<floatT> & _energDens;

              // vectors for the calculated angle and integrated energy density
              std::vector<floatT> _AngleSec1;
              std::vector<floatT> _EDensSec1;
              std::vector<floatT> _AngleSec2;
              std::vector<floatT> _EDensSec2;
              std::vector<floatT> _AngleSec3;
              std::vector<floatT> _EDensSec3;
              std::vector<floatT> _AngleSec4;
              std::vector<floatT> _EDensSec4;

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

       public:
              // constructor
              IntegratedEnergyDensity(Site startSite, EnergyDensity<floatT> & newEnergDens):
                                                        _RealX(startSite.x()),
                                                        _RealY(startSite.y()),
                                                        _JetStartX((int) ((startSite.x() + 15. ) * 100.) ),
                                                        _JetStartY((int) ((startSite.y() + 15. ) * 100.) ),
                                                        _AngleStep((PI / 2.0) / ((floatT) Steps)),
                                                        _energDens(newEnergDens),
                                                        _AngleSec1(Steps),
                                                        _EDensSec1(Steps),
                                                        _AngleSec2(Steps),
                                                        _EDensSec2(Steps),
                                                        _AngleSec3(Steps),
                                                        _EDensSec3(Steps),
                                                        _AngleSec3(Steps),
                                                        _EDensSec3(Steps),
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


              // bilinear interpolation function between the grid points in the energy density grid
              floatT energyDensityInterpolation(floatT fP11, floatT fP12, floatT fP21, floatT fP22,
                                          floatT x1, floatT x2, floatT x, floatT y1, floatT y2, floatT y){
              	floatT x2x, x2x1, xx1, y2y, y2y1, yy1, interpolationValue;
              	x2x = x2 - x;
              	x2x1 = x2 - x1;
              	xx1 = x - x1;
              	y2y = y2 - y;
              	y2y1 = y2 - y1;
              	yy1 = y - y1;
              	return interpolationValue = (1.0 / (x2x1 * y2y1)) * (fP11 * x2x * y2y + fP21 * xx1 * y2y
                                                 + fP12 * x2x * yy1 + fP22 * xx1 * yy1);
              }


              // calculation for the right sector (sector 1) from 7/4 pi to 1/4 pi
              void sector1() {
                     // PhiOne is the angle with respect to the origin in real space
                     floatT PhiOne = 0.0;
                     // alpha is the angle seen from the start point of the integration
                     floatT Alpha1 = -PI / 4.0;
                     for(int i = 0; i < Steps; i++) {
                            // calculate the slope from alpha
                     	floatT m = std::tan(Alpha1);
                     	floatT x1, x2, y, y1, y2;
                            // point-slope form of equation for straight line y = m(x-x1)+y1
                            // calculate for each x value of the grid a value for y
                     	for(int x = _JetStartX; x <= 3000; x++) {
                     		y = m * ((floatT) x - (floatT) _JetStartX) + (floatT) _JetStartY;
                                   // check wheather the point lies on the grid or not (needed for
                                   // integration which does not start at the origin (real (0.0,0.0))
                     		if(y <= 3000.0 && y >= 0.0) {
                                          // find the four gridpoints {(x1,y1),(x2,y1),(x1,y2),(x2,y2)}
                                          // for the interpolation around a point (x,y)
                     			x1 = std::floor((floatT) x);
                     			x2 = std::floor((floatT) x + 1.0);
                     			y1 = std::floor(y);
                     			y2 = std::floor(y + 1.0);
                                          // define the points where the energy density grid is read out
                     			Site siteP1((int) x1, (int) y1);
                                          Site siteP2((int) x1, (int) y2);
                                          Site siteP3((int) x2, (int) y1);
                                          Site siteP4((int) x2, (int) y2);

                     			floatT Interpol = energyDensityInterpolation(
                                                 _energDens.getSmearedEnergyDensData() -> getSite(siteP1),
                                                 _energDens.getSmearedEnergyDensData() -> getSite(siteP2),
                                                 _energDens.getSmearedEnergyDensData() -> getSite(siteP3),
                                                 _energDens.getSmearedEnergyDensData() -> getSite(siteP4),
                                                 x1, x2, (floatT) x, y1, y2, y);
                                          // calculate the distance from the jet origin
                     			floatT r = std::hypot((((floatT) x - (floatT) _JetStartX) / 100.0),
                                                        ((y - (floatT) _JetStartY) / 100.0));
                                          // store point at line and the interpolated value at that point
                     			_RadiusOne.at(x - _JetStartX) = r;
                     			_EDensOne.at(x - _JetStartX) = Interpol;

                     		}
                     		else {
                                          // calculate the distance from the jet origin
                     			floatT r = std::hypot((((floatT) x - (floatT) _JetStartX) / 100.0),
                                                 ((y - (floatT) _JetStartY) / 100.0));
                                          // store point at line and the interpolated value at that point
                                          _RadiusOne.at(x - _JetStartX) = r;
                     			_EDensOne.at(x - _JetStartX) = 0.0;
                     		}

                     	}
                            // perform a trapezoidal integration along the line
                     	floatT Integral = 0.0;
                     	for(int j = 0; j <= 3000 - _JetStartX; j++)
                     	{
                     		Integral += 0.5*(_RadiusOne.at(j+1) - _RadiusOne.at(j))
                                          * (_EDensOne.at(j) + _EDensOne.at(j+1));
                     	}

                     // calculate the angle from the center of mass (real (0.0,0.0))
                     y = m * (3000.0 - (floatT) _JetStartX) + (floatT) _JetStartY;
                     PhiOne = std::atan((y - 1500.0) / (1500.0));
                     if(PhiOne < 0.0) PhiOne += 2.0 * PI;

                     // store the values for the angle and the integral vectors
                     _AngleSec1.at(i) = PhiOne;
                     _EDensSec1.at(i) = Integral;

                     // Calculation of ECCENTRICITY

                     Alpha1 += _AngleStep;
                     }
              }

              // calculation for the left sector (sector 2) from 3/4 pi to 5/4 pi
              void sector2() {
                     // PhiTwo is the angle with respect to the origin in real space
                     floatT PhiTwo = 0.0;
                     // alpha is the angle seen from the start point of the integration
                     floatT Alpha2 = (3.0 * PI) / 4.0;
                     for(int i = 0; i < Steps; i++) {
                            // calculate the slope from alpha
                     	floatT m = std::tan(Alpha2);
                     	floatT x1, x2, y, y1, y2;
                            // point-slope form of equation for straight line y = m(x-x1)+y1
                            // calculate for each x value of the grid a value for y
                     	for(int x = _JetStartX; x >= 0; x--) {
                     		floatT y = m*((floatT) x - (floatT) _JetStartX) + (floatT) _JetStartY;
                                   // check wheather the point lies on the grid or not (needed for
                                   // integration which does not start at the origin (real (0.0,0.0))
                     		if(y <= 3000.0 && y >= 0.0) {
                                          // find the four gridpoints {(x1,y1),(x2,y1),(x1,y2),(x2,y2)}
                                          // for the interpolation around a point (x,y)
                     			x1 = std::floor((floatT) x);
                     			x2 = std::floor((floatT) x + 1.0);
                     			y1 = std::floor(y);
                     			y2 = std::floor(y + 1.0);
                                          // define the points where the energy density grid is read out
                                          Site siteP1((int) x1, (int) y1);
                                          Site siteP2((int) x1, (int) y2);
                                          Site siteP3((int) x2, (int) y1);
                                          Site siteP4((int) x2, (int) y2);

                                          floatT Interpol = energyDensityInterpolation(
                                                 _energDens.getSmearedEnergyDensData() -> getSite(siteP1),
                                                 _energDens.getSmearedEnergyDensData() -> getSite(siteP2),
                                                 _energDens.getSmearedEnergyDensData() -> getSite(siteP3),
                                                 _energDens.getSmearedEnergyDensData() -> getSite(siteP4),
                                                 x1, x2, (floatT) x, y1, y2, y);
                                          // calculate the distance from the jet origin
                     			floatT r = std::hypot((((floatT) x - (floatT) _JetStartX) / 100.0),
                                                        ((y - (floatT) _JetStartY) / 100.0));
                                          // store point at line and the interpolated value at that point
                     			_RadiusTwo.at(_JetStartX - x) = r;
                     			_EDensTwo.at(_JetStartX - x) = Interpol;
                     		}
                     		else {
                                          // calculate the distance from the jet origin
                                          floatT r = std::hypot((((floatT) x - (floatT) _JetStartX) / 100.0),
                                                 ((y - (floatT) _JetStartY) / 100.0));
                                          // store point at line and the interpolated value at that point
                                          _RadiusTwo.at(_JetStartX - x) = r;
                                          _EDensTwo.at(_JetStartX - x) = 0.0;
                     		}
                     	}

                            // perform a trapezoidal integration along the line
                            floatT Integral = 0.0;
                            for(int j = 0; j <= _JetStartX; j++)
                            {
                                   Integral += 0.5*(_RadiusTwo.at(j + 1) - _RadiusTwo.at(j))
                                          * (_EDensTwo.at(j) + _EDensTwo.at(j + 1));
                            }

                     // calculate the angle from the center of mass (real (0.0,0.0))
                     y = m * (0.0 - (floatT) _JetStartX) + (floatT) _JetStartY;
                     PhiTwo = atan((1500.0 - y) / (1500.0)) + PI;
                     if(PhiTwo < 0.0) PhiTwo += 2.0 * PI;

                     // store the values for the angle and the integral in an array
                     _AngleSec2.at(i) = PhiTwo;
                     _EDensSec2.at(i) = Integral;

                     // Calculation of ECCENTRICITY

                     Alpha2 += _AngleStep;
                     }
              }

              // calculation for the upper sector (sector3) from 1/4 pi to 3/4 pi
              void sector3() {
                     // the roles of x and y must be changed because the slope/tan diverges at pi/2
                     // PhiThree is the angle with respect to the origin in real space
                     floatT PhiThree = 0.0;

                     // alpha is the angle seen from the start point of the integration
                     floatT Alpha3 = -PI / 4.0;
                     for(int i = 0; i < Steps; i++) {
                            // calculate the slope from alpha
                     	floatT m = std::tan(Alpha3);
                     	floatT x1, x2, x, y1, y2;
                            // point-slope form of equation for straight line x = m(y-y1)+x1
                            // calculate for each y value of the grid a value for x
                     	for(int y = _JetStartY; y <= 3000; y++) {
                     		floatT x = -m * ((floatT) y - (floatT) _JetStartY) + (floatT) _JetStartX;
                                   // check wheather the point lies on the grid or not (needed for
                                   // integration which does not start at the origin (real (0.0,0.0))
                     		if(x <= 3000.0 && x >= 0.0) {
                                          // find the four gridpoints {(x1,y1),(x2,y1),(x1,y2),(x2,y2)}
                                          // for the interpolation around a point (x,y)
                     			y1 = std::floor((floatT) y);
                     			y2 = std::floor((floatT) y + 1.0);
                     			x1 = std::floor(x);
                     			x2 = std::floor(x + 1.0);
                                          // define the points where the energy density grid is read out
                                          Site siteP1((int) x1, (int) y1);
                                          Site siteP2((int) x1, (int) y2);
                                          Site siteP3((int) x2, (int) y1);
                                          Site siteP4((int) x2, (int) y2);

                                          floatT Interpol = energyDensityInterpolation(
                                                 _energDens.getSmearedEnergyDensData() -> getSite(siteP1),
                                                 _energDens.getSmearedEnergyDensData() -> getSite(siteP2),
                                                 _energDens.getSmearedEnergyDensData() -> getSite(siteP3),
                                                 _energDens.getSmearedEnergyDensData() -> getSite(siteP4),
                                                 x1, x2, x, y1, y2, (floatT) y);
                                          // calculate the distance from the jet origin
                            		floatT r = std::hypot((((floatT) x - (floatT) _JetStartX) / 100.0),
                                                        ((y - (floatT) _JetStartY) / 100.0));
                                          // store point at line and the interpolated value at that point
                     			_RadiusThree.at(y - _JetStartY) = r;
                     			_EDensThree.at(y - _JetStartY) = Interpol;
                     		}
                     		else {
                                          // calculate the distance from the jet origin
                            		floatT r = std::hypot((((floatT) x - (floatT) _JetStartX) / 100.0),
                                                        ((y - (floatT) _JetStartY) / 100.0));
                                          // store point at line and the interpolated value at that point
                     			_RadiusThree.at(y - _JetStartY) = r;
                     			_EDensThree.at(y - _JetStartY) = 0.0;
                     		}
                     	}
                            // perform a trapezoidal integration along the line
                     	floatT Integral = 0.0;
                     	for(int j = 0; j <= 3000 - _JetStartY; j++)
                     	{
                     		Integral += 0.5 * (_RadiusThree.at(j + 1) - _RadiusThree.at(j))
                                          * (_EDensThree.at(j) + _EDensThree.at(j + 1));
                     	}

                     // calculate the angle from the center of mass (real (0.0,0.0))
                     x = m * (3000.0 - (floatT) _JetStartY) + (floatT) _JetStartX;
                     PhiThree = std::atan((x - 1500.0) / (1500.0)) + PI / 2.0;
                     if(PhiThree < 0.0) PhiThree += 2.0 * PI;

                     // store the values for the angle and the integral in an array
                     _AngleSec3.at(i) = PhiThree;
                     _EDensSec3.at(i) = Integral;

                     // Calculation of ECCENTRICITY

                     Alpha3 += _AngleStep;
                     }
              }

              // calculation for the lower sector from 5/4 pi to 7/4 pi
              void sector4() {
                     // the roles of x and y must be changed because the slope/tan diverges at pi/2
                     // PhiThree is the angle with respect to the origin in real space
                     floatT PhiFour = 0.0;

                     // alpha is the angle seen from the start point of the integration
                     floatT Alpha4 = -PI / 4.0;
                     for(int i = 0; i < Steps; i++)
                     {
                            // calculate the slope from alpha
                     	floatT m = std::tan(Alpha4);
                     	floatT x1, x2, x, y1, y2;
                     //point-slope form of equation for straight line x = m(y-y1)+x1
                     //calculate for each y value of the grid a value for x
                     	for(int y = _JetStartY; y >= 0; y--)
                     	{
                     		floatT x = -m * ((floatT) y - (floatT) _JetStartY) + (floatT) _JetStartX;
                                   // check wheather the point lies on the grid or not (needed for
                                   // integration which does not start at the origin (real (0.0,0.0))
                     		if(x <= 3000.0 && x >= 0.0) {
                                          // find the four gridpoints {(x1,y1),(x2,y1),(x1,y2),(x2,y2)}
                                          // for the interpolation around a point (x,y)
                     			y1 = std::floor((floatT) y);
                     			y2 = std::floor((floatT) y + 1.0);
                     			x1 = std::floor(x);
                     			x2 = std::floor(x + 1.0);
                                          // define the points where the energy density grid is read out
                                          Site siteP1((int) x1, (int) y1);
                                          Site siteP2((int) x1, (int) y2);
                                          Site siteP3((int) x2, (int) y1);
                                          Site siteP4((int) x2, (int) y2);

                                          floatT Interpol = energyDensityInterpolation(
                                                 _energDens.getSmearedEnergyDensData() -> getSite(siteP1),
                                                 _energDens.getSmearedEnergyDensData() -> getSite(siteP2),
                                                 _energDens.getSmearedEnergyDensData() -> getSite(siteP3),
                                                 _energDens.getSmearedEnergyDensData() -> getSite(siteP4),
                                                 x1, x2, x, y1, y2, (floatT) y);
                                          // calculate the distance from the jet origin
                            		floatT r = std::hypot((((floatT) x - (floatT) _JetStartX) / 100.0),
                                                        ((y - (floatT) _JetStartY) / 100.0));
                                          // store point at line and the interpolated value at that point

                     			_RadiusFour.at(_JetStartY - y) = r;
                     			_EDensFour.at(_JetStartY - y) = Interpol;
                     		}
                     		else {
                                          // calculate the distance from the jet origin
                            		floatT r = std::hypot((((floatT) x - (floatT) _JetStartX) / 100.0),
                                                        ((y - (floatT) _JetStartY) / 100.0));
                                          // store point at line and the interpolated value at that point

                     			_RadiusFour.at(_JetStartY - y) = r;
                     			_EDensFour.at(_JetStartY - y) = 0.0;
                     		}
                     	}
                            //perform a trapezoidal integration along the line
                     	floatT Integral = 0.0;
                     	for(int j = 0; j <= _JetStartY; j++)
                     	{
                     		Integral += 0.5 * (_RadiusFour.at(j + 1) - _RadiusFour.at(j))
                                          * (_EDensFour.at(j) + _EDensFour.at(j + 1));
                     	}

                     //calculate the angle from the center of mass (real (0.0,0.0))
                     x = m * (0.0 - (floatT) _JetStartY) + (floatT) _JetStartX;
                     PhiFour = std::atan((1500.0 - x) / (1500.0)) + (6.0 * PI) / 4.0;
                     if(PhiFour < 0.0) PhiFour += 2.0 * PI;

                     // store the values for the angle and the integral in an array
                     _AngleSec4.at(i) = PhiFour;
                     _EDensSec4.at(i) = Integral;

                     // Calculation of ECCENTRICITY

                     Alpha4 += _AngleStep;
                     }
              }

              floatT inline integralNormalization() {
                     // normalization of the integral by summing the value of the integrals
                     // for all angles phi of one event
                     floatT IntegralNormalization = 0.0;
                     for(int i = 0; i < Steps; i++)
                     {
                     	IntegralNormalization += (_EDensSec1.at(i) + _EDensSec2.at(i)
                                                 + _EDensSec3.at(i) + _EDensSec4.at(i));
                     }
                     return IntegralNormalization;
              }

              void averagedIntegral() {
                     floatT IntegralNormalization = integralNormalization();
                     // divide each integral by the normalization and the number of events
                     for(int j = 0; j < Steps; j++)
                     {
                     	_AverageEDens1.at(j) += _EDensSec1.at(j) / IntegralNormalization;
                     	_AverageEDens2.at(j) += _EDensSec2.at(j) / IntegralNormalization;
                     	_AverageEDens3.at(j) += _EDensSec3.at(j) / IntegralNormalization;
                     	_AverageEDens4.at(j) += _EDensSec4.at(j) / IntegralNormalization;
                     }
              }

              friend class Eccentricity;
              friend class FlowCoefficients;
};


template<class floatT>
class Eccentricity{
       private:
              // vectors for the eccentricity values of all sectors
              std::vector<std::complex<floatT>> _EccCounter1;
              std::vector<std::complex<floatT>> _EccDenom1;

              std::vector<std::complex<floatT>> _EccCounter2;
              std::vector<std::complex<floatT>> _EccDenom2;

              std::vector<std::complex<floatT>> _EccCounter3;
              std::vector<std::complex<floatT>> _EccDenom3;

              std::vector<std::complex<floatT>> _EccCounter4;
              std::vector<std::complex<floatT>> _EccDenom4;

              std::vector<floatT> _EventEccentricity;

              IntegratedEnergyDensity<floatT> & _intEnergDens;

       public:
              // constructor
              Eccentricity(IntegratedEnergyDensity<floatT> & newIntEnergDens):
                            _intEnergDens(newIntEnergDens),
                            _EccCounter1(4),
                            _EccDenom1(4),
                            _EccCounter2(4),
                            _EccDenom2(4),
                            _EccCounter3(4),
                            _EccDenom3(4),
                            _EccCounter4(4),
                            _EccDenom4(4),
                            _EventEccentricity(4) {}

              // calculation of the eccentricity for sector 1
              void eccentricitySector1() {
                     //define the imaginary I
                     const std::complex<floatT> I(0.0,1.0);
                     // perform a trapezoidal integration to calculate the eccentricity e2 ... e5
                     for (int i = 0; i < Steps; i++) {
                            for(int n =  2; n <= 5; n++) {
                                   for(int j = 0; j < 3000 - _intEnergDens._JetStartX; j++) {
                                          _EccCounter1.at(n - 2) += (_intEnergDens._RadiusOne.at(j + 1)
                                                        - _intEnergDens._RadiusOne.at(j)) * 0.5
                                                        * (_intEnergDens._EDensOne.at(j)
                                                        * std::pow(_intEnergDens._RadiusOne.at(j),(n + 1))
                                                        + _intEnergDens._EDensOne.at(j + 1)
                                                        * std::pow(_intEnergDens._RadiusOne.at(j + 1),(n + 1)))
                                                        * exp(I * ((floatT) n * _intEnergDens._AngleSec1.at(i)));

                                          _EccDenom1.at(n - 2) += (_intEnergDens._RadiusOne.at(j + 1)
                                                        - _intEnergDens._RadiusOne.at(j)) * 0.5
                                                        * (_intEnergDens._EDensOne.at(j)
                                                        * std::pow(_intEnergDens._RadiusOne.at(j),(n + 1))
                                                        + _intEnergDens._EDensOne.at(j + 1)
                                                        * std::pow(_intEnergDens._RadiusOne.at(j + 1),(n + 1)));
                                   }
                            }
                     }
              }

              void eccentricitySector2() {
                     //define the imaginary I
                     const std::complex<floatT> I(0.0,1.0);
                     // perform a trapezoidal integration to calculate the eccentricity e2 ... e5
                     for (int i = 0; i < Steps; i++) {
                            for(int n =  2; n <= 5; n++) {
                                   for(int j = 0; j < 3000 - _intEnergDens._JetStartX; j++){
                                          _EccCounter2.at(n - 2) += (_intEnergDens._RadiusTwo.at(j + 1)
                                                 - _intEnergDens._RadiusTwo.at(j)) * 0.5
                                                 * (_intEnergDens._EDensTwo.at(j)
                                                 * std::pow(_intEnergDens._RadiusTwo.at(j),(n + 1))
                                                 + _intEnergDens._EDensTwo.at(j + 1)
                                                 * std::pow(_intEnergDens._RadiusTwo.at(j + 1),(n + 1)))
                                                 * exp(I * ((floatT) n * _intEnergDens._AngleSec2.at(i)));

                                          _EccDenom2.at(n - 2) += (_intEnergDens._RadiusTwo.at(j + 1)
                                                 - _intEnergDens._RadiusTwo.at(j)) * 0.5
                                                 * (_intEnergDens._EDensTwo.at(j)
                                                 * std::pow(_intEnergDens._RadiusTwo.at(j),(n + 1))
                                                 +_intEnergDens._EDensTwo.at(j + 1)
                                                 * std::pow(_intEnergDens._RadiusTwo.at(j + 1),(n + 1)));
                                   }
                            }
                     }
              }

              void eccentricitySector3() {
                     //define the imaginary I
                     const std::complex<floatT> I(0.0,1.0);
                     // perform a trapezoidal integration to calculate the eccentricity e2 ... e5
                     for (int i = 0; i < Steps; i++) {
                            for(int n =  2; n <= 5; n++) {
                                   for(int j = 0; j < 3000 - _intEnergDens._JetStartX; j++) {
                                          _EccCounter3.at(n - 2) += (_intEnergDens._RadiusThree.at(j + 1)
                                                        - _intEnergDens._RadiusThree.at(j)) * 0.5
                                                        * (_intEnergDens._EDensThree.at(j)
                                                        * std::pow(_intEnergDens._RadiusThree.at(j),(n + 1))
                                                        + _intEnergDens._EDensThree.at(j + 1)
                                                        * std::pow(_intEnergDens._RadiusThree.at(j + 1),(n + 1)))
                                                        * exp(I * ((floatT) n * _intEnergDens._AngleSec3.at(i)));

                                          _EccDenom3.at(n - 2) += (_intEnergDens._RadiusThree.at(j + 1)
                                                        - _intEnergDens._RadiusThree.at(j)) * 0.5
                                                        * (_intEnergDens._EDensThree.at(j)
                                                        * std::pow(_intEnergDens._RadiusThree.at(j),(n + 1))
                                                        + _intEnergDens._EDensThree.at(j + 1)
                                                        * std::pow(_intEnergDens._RadiusThree.at(j + 1),(n + 1)));
                                   }
                            }
                     }

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


                     //calculate the eccentricity e2 ... e5 for one event and store the value in the array
                     for(int n = 0; n < 4; n++){
              	             _EventEccentricity.at(n) = std::abs(-(_EccCounter1.at(n) + _EccCounter2.at(n)
                                          + _EccCounter3.at(n) + _EccCounter4.at(n))
                                          / (_EccDenom1.at(n) + _EccDenom2.at(n)
                                          + _EccDenom3.at(n) + _EccDenom4.at(n)));
                     }
              }

};


template<class floatT>
class FlowCoefficients{
       private:
              std::vector<floatT> _MergedAngles;
              std::vector<floatT> _MergedIntegrals;

              std::vector<floatT> _Fourier;

              IntegratedEnergyDensity<floatT> & _intEnergDens;

       public:
              // constructor
              FlowCoefficients(EnergyDensity<floatT> & newEnergDens):
                            _intEnergDens(newEnergDens),
                            _MergedAngles(4 * Steps),
                            _MergedIntegrals(4 * Steps),
                            _Fourier(6) {}

              // merged vectors
              void mergeVectors() {
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
              void flowCoefficients() {
                     floatT a0 = 0.0;
                     for(int i=0; i < 4 * Steps - 1; i++) {
                     	a0 += (((_MergedIntegrals.at(i) + _MergedIntegrals.at(i + 1)) / 2.0)
                                   * _intEnergDens._AngleStep) / PI;
                     }

                     floatT aFourier = 0.0;
                     for(int aCount=1; aCount <= 3; aCount++) {
                     	aFourier = 0.0;
                     	for(int i=0; i < 4 * Steps - 1; i++) {
                     		aFourier += ((((_MergedIntegrals.at(i) + _MergedIntegrals.at(i + 1)) / 2.0)
                                                 * std::cos((floatT) aCount
                                                               * ((_MergedAngles.at(i)
                                                               + _MergedAngles.at(i + 1)) / 2.0)))
                                                 * _intEnergDens._AngleStep) / PI;
                     	}
                     	_Fourier.at(aCount - 1) = aFourier / (a0 / 2.0);
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
                     	_Fourier.at(bCount + 2) = bFourier / (a0 / 2.0);
                     	bFourier = 0.0;
                     }
                     a0 = 0.0;
              }
};


int main(int argc, char const *argv[]) {
       std::clock_t start;

       start = std::clock();


       // runThroughGridTest(3);

       computeTest();

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
       std::cout << "now at (x,y): (" << site.x() << "," << site.y() << ")" << '\n';
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

       FileWriter<PREC> file(NEvents,NNucleonsCore, filename);

       file.readFile(rawDataGrid);

       std::cout << "Compute energy density" << '\n';

       energDens.energyDensity(rawDataGrid);

       std::cout << "Smear energy density" << '\n';

       energDens.smearedEnergyDensity();

       file.writeFileGrid(energDens.getSmearedEnergyDensData());

       std::cout << "Compute integrations in all directions" << '\n';

       IntegratedEnergyDensity<PREC> intEnergDens(Site(0, 0), energDens);

       intEnergDens.sector1();
       intEnergDens.sector2();
       intEnergDens.sector3();
       intEnergDens.sector4();

       intEnergDens.averagedIntegral();
}
