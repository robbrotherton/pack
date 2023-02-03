#include <Rcpp.h>
#include <math.h>
#include <random>
#include <algorithm>
using namespace Rcpp;
using namespace std;

// this returns TRUE if point is within boundary
bool polygonContainsPoint(NumericVector polyXPoints,
                          NumericVector polyYPoints,
                          double testX,
                          double testY);

// declare vecmin function
double vecmin(NumericVector x);
double vecmax(NumericVector x);

//' Pack circles
//'
//' This function returns a dataframe of x and y coordinates and radii
//' of circles packed into an arbitrary polygon container.
//'
//' @param polygon A dataframe containing the x and y coordinates of vertices
//'   a polygon into which circles will be packed
//' @param radii A numeric vector giving the radii of circles which will be
//'   packed into the container specified in the polygon argument. This also
//'   determines the maximum number of circles which the algorithm will attempt
//'   to place.
//' @param existing_circles A dataframe of x and y coordinates and radii of any
//'   existing circles to be avoided when packing the new ones. This allows
//'   multiple iterations, for example to pack a shape within another shape.
//' @param max_attempts A scalar integer. How many times the algorithm will try
//'   to place a circle in a new spot before giving up.
//' @param seed A scalar integer. Used to seed the random number generator.
//' @param neat_edges Logical. If true, no part of a circle can extend beyond
//'   the container polygons boundary. If false, circles can extend beyond the
//'   boundary as long as any part of the circle remains within the container.
//' @export
// [[Rcpp::export]]
DataFrame pack_circles(DataFrame polygon,
                       NumericVector radii,
                       DataFrame existing_circles,
                       int max_attempts = 2000,
                       int seed = 1,
                       bool neat_edges = true) {

  NumericVector x_existing = existing_circles["x"];
  NumericVector y_existing = existing_circles["y"];
  NumericVector r_existing = existing_circles["r"];
  int n_existing = x_existing.size();

  int n = 0;
  int max_circles = radii.size();

  int attempt = 0;
  bool overlap;
  double new_x;
  double new_y;
  double new_r;
  double dist;
  NumericVector polygonXPoints = polygon["x"];
  NumericVector polygonYPoints = polygon["y"];
  double minX = vecmin(polygonXPoints);
  double maxX = vecmax(polygonXPoints);
  double minY = vecmin(polygonYPoints);
  double maxY = vecmax(polygonYPoints);
  // NumericVector::iterator minY = std::min_element(polygonYPoints.begin(), polygonYPoints.end());
  // NumericVector::iterator maxX = std::max_element(polygonXPoints.begin(), polygonXPoints.end());
  // NumericVector::iterator maxY = std::max_element(polygonYPoints.begin(), polygonYPoints.end());
  // int numVerts = polygonXPoints.length();
  // int j = numVerts - 1;
  NumericVector x_out(max_circles);
  NumericVector y_out(max_circles);
  NumericVector r_out(max_circles);
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> rand_x(minX,maxX);
  std::uniform_real_distribution<double> rand_y(minY,maxY);

  while(attempt <= max_attempts) {

    overlap = FALSE;
    bool within_perim = FALSE;


    // first, make sure the new circle is placed within the polygon boundary
    while(within_perim==FALSE) {

      ++attempt;

      // if(attempt > max_attempts) break;
      new_x = rand_x(generator);
      new_y = rand_y(generator);
      new_r = radii[n];

      if(neat_edges) {

        for(int k = 0; k < 16; k++) {
          double theta = M_PI * (k / 8.0);
          double pointX = new_x + new_r * cos(theta);
          double pointY = new_y + new_r * sin(theta);
          if(!polygonContainsPoint(polygonXPoints,
                                   polygonYPoints,
                                   pointX,
                                   pointY)) {
            break;
          }
          if(k==15) within_perim=TRUE;
        }
      }

      if(!neat_edges) {
        within_perim = polygonContainsPoint(polygonXPoints,
                                            polygonYPoints,
                                            new_x,
                                            new_y);
      }
    }

    // now check the values against existing circles
    // returns TRUE if circles overlap

    for(int k = 0; k < n_existing; ++k) {

      dist = sqrt(pow(x_existing[k]-new_x, 2) + pow(y_existing[k]-new_y, 2));

      overlap = dist < (new_r + r_existing[k]);
      // if there's overlap, stop checking
      if(overlap) break;

    }

    if(!overlap) {
      for(int i = 0; i < n; ++i) {

        dist = sqrt(pow(x_out[i]-new_x, 2) + pow(y_out[i]-new_y, 2));

        overlap = dist < (new_r + r_out[i]);
        // if there's overlap, stop checking
        if(overlap) break;

      }}

    // if there's no overlap after checking every row...
    if(!overlap) {

      x_out[n] = new_x;
      y_out[n] = new_y;
      r_out[n] = new_r;

      ++n;
      attempt = 0;

    }

    if(n >= max_circles) break; // all the circles fit

  }

  Rcout << "made " << n;
  x_out.erase(std::remove(x_out.begin(), x_out.end(), 0), x_out.end());
  y_out.erase(std::remove(y_out.begin(), y_out.end(), 0), y_out.end());
  r_out.erase(std::remove(r_out.begin(), r_out.end(), 0), r_out.end());

  return DataFrame::create(_["x"]= x_out, _["y"]= y_out, _["r"]= r_out);

}



// [[Rcpp::export]]
SEXP pack_polygons(DataFrame polygon,
                   NumericVector radii,
                   int sides,
                   int max_attempts = 2000,
                   int seed = 1,
                   bool neat_edges = true) {

  // NumericVector x_existing = existing_circles["x"];
  // NumericVector y_existing = existing_circles["y"];
  // NumericVector r_existing = existing_circles["r"];
  // int n_existing = x_existing.size();

  int n = 0;
  int p = sides*(9-sides);
  int max_circles = radii.size();
  // int sides = 5;
  int attempt = 0;
  bool overlap;
  double new_x;
  double new_y;
  double new_r;
  // double theta;
  // double denominator;
  NumericVector verts_x(p+2);
  NumericVector verts_y(p+2);
  // NumericVector verts_x;
  // NumericVector verts_y;
  // double dist;
  NumericVector polygonXPoints = polygon["x"];
  NumericVector polygonYPoints = polygon["y"];
  double minX = vecmin(polygonXPoints);
  double maxX = vecmax(polygonXPoints);
  double minY = vecmin(polygonYPoints);
  double maxY = vecmax(polygonYPoints);
  NumericVector x_out(max_circles);
  NumericVector y_out(max_circles);
  NumericVector r_out(max_circles);
  List verts_x_out(max_circles);
  List verts_y_out(max_circles);
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> rand_x(minX,maxX);
  std::uniform_real_distribution<double> rand_y(minY,maxY);


  while(attempt <= max_attempts) {

    overlap = FALSE;
    bool within_perim = FALSE;


    // first, make sure the new circle is placed within the polygon boundary
    while(within_perim==FALSE) {

      ++attempt;

      // instead of coming up with center and radius, need to define the
      // entire polygon -- every vertex
      // int sides = 4;
      new_r = radii[n];
      new_x = rand_x(generator);
      new_y = rand_y(generator);
      verts_x = NumericVector(p+2);
      verts_y = NumericVector(p+2);
      // verts_x = Rcpp::NumericVector::create(0,0,0,0,0);
      // verts_y = Rcpp::NumericVector::create(0,0,0,0,0);
      // RawVector verts_y(sides);

      for(int v = 0; v <= p+1; v++) {
        // double theta = PI * v / sides * 2 + PI/4;
        double theta = M_PI*2*v/p;
        double numerator = cos(M_PI/sides);
        double denominator = cos(2.0/sides * asin(cos(sides/2.0*theta)));
        double r = numerator/denominator;
        // Rcout << theta << "//" << numerator << "/" << denominator << "\n";

        verts_x[v] = new_x + r * new_r * cos(theta);
        verts_y[v] = new_y + r * new_r * sin(theta);
      }

      if(neat_edges) {

        for(int k = 0; k < verts_x.size(); k++) {
          if(!polygonContainsPoint(polygonXPoints,
                                   polygonYPoints,
                                   verts_x[k],
                                          verts_y[k])) {
            break;
          }
          if(k==verts_x.size()-1) within_perim=TRUE;
        }
      }

      if(!neat_edges) {
        within_perim = polygonContainsPoint(polygonXPoints,
                                            polygonYPoints,
                                            new_x,
                                            new_y);
      }
    }

    // now check the values against existing circles
    // returns TRUE if circles overlap
    //
    //     for(int k = 0; k < n_existing; ++k) {
    //
    //       dist = sqrt(pow(x_existing[k]-new_x, 2) + pow(y_existing[k]-new_y, 2));
    //
    //       overlap = dist < (new_r + r_existing[k]);
    //
    //
    //       // if there's overlap, stop checking
    //       if(overlap) break;

    //     }

    // if there's no overlap with existing circles, check the new ones...,
    // if(!overlap) {
    for(int i = 0; i < n; ++i) {

      // NumericVector theta

      // need to interpolate more points between vertices...
      for(int v = 0; v < verts_x.size(); ++v) {
        overlap = polygonContainsPoint(verts_x_out[i],
                                       verts_y_out[i],
                                                  verts_x[v],
                                                         verts_y[v]);
        if(overlap) break;
      }
      //
      //         dist = sqrt(pow(x_out[i]-new_x, 2) + pow(y_out[i]-new_y, 2));
      //
      //         overlap = dist < (new_r + r_out[i]);
      // if there's overlap, stop checking
      if(overlap) break;

    }
    // }

    // if there's no overlap after checking every row...
    if(!overlap) {

      x_out[n] = new_x;
      y_out[n] = new_y;
      r_out[n] = new_r;
      verts_x_out[n] = verts_x;
      verts_y_out[n] = verts_y;

      ++n;
      attempt = 0;

    }

    if(n >= max_circles) break; // all the circles fit

  }

  Rcout << "made " << n;
  x_out.erase(std::remove(x_out.begin(), x_out.end(), 0), x_out.end());
  y_out.erase(std::remove(y_out.begin(), y_out.end(), 0), y_out.end());
  r_out.erase(std::remove(r_out.begin(), r_out.end(), 0), r_out.end());
  verts_x_out.erase(std::remove(verts_x_out.begin(), verts_x_out.end(), R_NilValue), verts_x_out.end());
  verts_y_out.erase(std::remove(verts_y_out.begin(), verts_y_out.end(), R_NilValue), verts_y_out.end());
  // return DataFrame::create(_["x"]= x_out, _["y"]= y_out, _["r"]= r_out);


  // Construct a list from the columns
  Rcpp::List ret = Rcpp::List::create(
    Rcpp::Named("x") = x_out,
    Rcpp::Named("y") = y_out,
    Rcpp::Named("r") = r_out,
    Rcpp::Named("x_verts") = verts_x_out,
    Rcpp::Named("y_verts") = verts_y_out);

  // Coerce to a data.frame
  ret.attr("class") = "data.frame";
  ret.attr("row.names") = Rcpp::seq(1, x_out.size());

  // Return the data.frame
  return ret;


}

// define the perimeter checking function here
bool polygonContainsPoint(NumericVector polyXPoints,
                          NumericVector polyYPoints,
                          double testX,
                          double testY) {
  int numVerts = polyXPoints.length();
  bool c = false;
  int j = numVerts - 1;
  for (int i = 0; i < numVerts; i++)
  {
    double deltaX = polyXPoints[j] - polyXPoints[i];
    double ySpread = testY - polyYPoints[i];
    double deltaY = polyYPoints[j] - polyYPoints[i];
    if (((polyYPoints[i] > testY) != (polyYPoints[j] > testY)) &&
        (testX < (((deltaX * ySpread) / deltaY) + polyXPoints[i])))
    {
      c = !c;
    }

    j = i;
  }
  return c;
}


// // function to get minimum element of vector
double vecmin(NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference
  return *it;
}


double vecmax(NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::max_element(x.begin(), x.end());
  // we want the value so dereference
  return *it;
}
