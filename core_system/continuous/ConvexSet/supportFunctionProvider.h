/**
 * Base class to represent a convex set whose support function can be computed. This is an ABSTRACT class.
 * Any convex set whose support function can be computed should be derived from this base class
 * @Rajarshi
 */

#ifndef __SUPPFUNC_PROV__
#define __SUPPFUNC_PROV__

#include <vector>
#include "core_system/math/matrix.h"
//#include "core_system/math/glpk_lp_solver/glpk_lp_solver.h"
#include "core_system/math/lp_solver/lp_solver.h"
#include <boost/shared_ptr.hpp>
#include "core_system/HybridAutomata/vartoindexmap.h"

class supportFunctionProvider : public var_to_index_map
{

public:
	typedef boost::shared_ptr<supportFunctionProvider> ptr;
	supportFunctionProvider() {
	}
	;
	virtual ~supportFunctionProvider() {
	}
	;
	/** Returns the dimension of the continuous set */
	virtual unsigned int getSystemDimension() const = 0;
	virtual bool getIsEmpty() const = 0;
	/**
	 * The compute support will be a function of the support function of the initial set
	 * and the input set.
	 */
	virtual double computeSupportFunction(std::vector<double> direction,
			lp_solver &lp) = 0;
	virtual double max_norm(int lp_solver_type_choosen,
			unsigned int dim_size) = 0;
};

#endif
