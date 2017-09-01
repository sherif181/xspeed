/*
 * Polytope.h
 *
 *  Created on: 30-Jun-2014
 *      Author: amit
 */

#ifndef POLYTOPE_H_
#define POLYTOPE_H_

#include <vector>
#include "core_system/math/matrix.h"
//#include "core_system/math/glpk_lp_solver/glpk_lp_solver.h"
#include "core_system/math/lp_solver/lp_solver.h"
#include "core_system/continuous/ConvexSet/supportFunctionProvider.h"
#include <boost/shared_ptr.hpp>
#include<set>
#include<utility>
#include <boost/tokenizer.hpp>

using namespace std;
/*
 * If a polytope is represented using intersections of halfspaces then it is of the form Ax<=b
 * where A is the coefficients Matrix of the variables 'x' and b is the columnVector
 *
 * coeffMatrix : 	All facets of the Polytope in Matrix form (i.e. the coefficients of the variables).
 * number_facets : Number of faces of the defined polytope.
 * dimension :  Number of variables of the system.
 * columnVector :	The values b for each facets.
 * InEqualitySign :	The in equalities sign of the bound values 'b'. Possible values are
 * 					0 :	for  Ax = b (b Equals to)
 * 					1 :	for  Ax <= b (b is Greater Than or Equals to)
 * 					2 :	for  Ax >= b (b is Less Than or Equals to)
 *
 * Also include function to compute Support Function for any given direction w.r.t. the defined polytope.
 * Also include function to return the Dimension of the defined polytope.
 */

class polytope: public supportFunctionProvider {

private:
	//glpk_lp_solver lp;	//Create only one lp at the time of creation of the polytope
	math::matrix<double> coeffMatrix;
	std::vector<double> columnVector;
	int InEqualitySign;
	unsigned int number_facets;
	unsigned int system_dimension;
	//bool lp_init;
	bool IsEmpty;
	bool IsUniverse;

public:
	typedef boost::shared_ptr<polytope> ptr;
	polytope();
	polytope(bool empty);
	polytope(math::matrix<double> coeffMatrix, std::vector<double> columnVector,
			int InEqualitySign);
	void setIsEmpty(bool empty);
	bool getIsEmpty() const;
	void setIsUniverse(bool universe);
	bool getIsUniverse();

	void setPolytope(math::matrix<double> coeffMatrix,
			std::vector<double> columnVector, int inEqualitySign);
	/*
	 * Adds one constraint to the existing polytope by adding the
	 * coefficient constraint with the bound value to the existing list.
	 */
	void setMoreConstraints(std::vector<double> coeff_constraint,
			double bound_value);

	/*
	 * Adds one or more constraints to the existing polytope by adding the
	 * coefficient_constraints with the bound_values to the existing list.
	 */
	void setMoreConstraints(math::matrix<double> coeff_constraints,
			std::vector<double> bound_values);

	// void set_Default_lp_init();
	// void set_lp_object(glpk_lp_solver* newObject);
	//const std::vector<std::vector<double> > getCoeffMatrix() const;
	math::matrix<double>& getCoeffMatrix();
	void setCoeffMatrix(const math::matrix<double> coeffMatrix);
	void setColumnVector(const std::vector<double> columnVector);

	int getInEqualitySign() const;
	void setInEqualitySign(int inEqualitySign);
	std::vector<double> getColumnVector();

	unsigned int getSystemDimension() const;
	void setSystemDimension(unsigned int systemDimension); //returns the number of variables of the polytopes.
	unsigned int getNumberFacets() const;
	void setNumberFacets(unsigned int numberFacets);

	double computeSupportFunction(std::vector<double> direction, lp_solver &lp);

	/*
	 * Returns Max norm of the polytope
	 */
	double max_norm(int lp_solver_type_choosen, unsigned int dim_for_Max_Norm);
	/*
	 * Returns True if polytope P1(the calling polytope object) and P2 intersects each other
	 *  i.e., True iff	P1 intersection P2 != empty set
	 */
	bool check_polytope_intersection(polytope::ptr P2,
			int lp_solver_type_choosen);
	/*
	 * Returns a new polytope after appending the constraints of P2
	 * which is an intersection-region
	 */
	const polytope::ptr GetPolytope_Intersection(polytope::ptr P2);

	/**
	 * Enumerate all vertices of the polytope between the two vectors
	 * given as arguments
	 */
	void enum_2dVert_restrict(std::vector<double> u, std::vector<double> v,
			int i, int j, std::set<std::pair<double, double> >&pts);

	/**
	 * enumerate all vertices of the polytope
	 * int i is the 1st projecting variable and
	 * 	   j is the second projecting variable
	 * 	   the value/index of i and j begins with 0 to n-1
	 */
	std::set<std::pair<double, double> > enumerate_2dVertices(int i, int j);


	/*
	 * Returns the list of vertices of the polytope in 2d with the given inputs as
	 * i and j where i and j are the 1st and 2nd projecting variables
	 */

	math::matrix<double> get_2dVertices(int dim1, int dim2);

	/**
	 * Computes the distance of a point from the polytope.
	 * If the point is inside the polytope, a 0 distance
	 * is returned. Otherwise, the distance is the sum of
	 * the all point to facet distances.
	 */
	double point_distance(std::vector<double> v);

	/*
	 * Prints the vertices of the polytope to a file, passed as parameter.
	 * The file could be called with any plotting utility.
	 *
	 */

	void print2file(std::string fname, unsigned int dim1, unsigned int dim2);

	/*
	 * debug function
	 */
	void print2files();

};


/**
 * Creates a pair of <loc_id, poly> from the user given bad state string
 */
void string_to_poly(const std::string& bad_state, std::pair<int, polytope::ptr>& f_set);


#endif /* POLYTOPE_H_ */
