
#include "preconditioner.h"
#include "solution.h"
#include "util.h"

#include <cmath>


/** ************************************************************************
 * Base constructor  for the Preconditioner class. 
 *
 *
 * @param size The number of grid points used in the approximation.
 * ************************************************************************ */
Preconditioner::Preconditioner(int number)
{
	setN(number);
	// allocate the vector with the lower diagonal matrix entries for
	// the Cholesky decomposition of the second order finite
	// difference operator for the same equation.
	vector = ArrayUtils<double>::twotensor(number,2);	
}

/** ************************************************************************
 *	Copy constructor  for the Preconditioner class. 
 *
 *	@param oldCopy The Preconditioner class member to make a copy of.
 * ************************************************************************ */
Preconditioner::Preconditioner(const Preconditioner& oldCopy)
{
	setN(oldCopy.getN());
	// Allocate the vector for the preconditioner, and copy it over.
	vector = ArrayUtils<double>::twotensor(getN(), Rank.COLNUMBER);
	int lupe;
	for(lupe=getN()-1;lupe>=0;lupe--)
		{
			for(int innerlupe= Rank.COLNUMBER-1; innerlupe>=0; innerlupe--)
			{
                vector[lupe][innerlupe] = oldCopy.getValue(lupe,innerlupe);
			}
			
		}
}

/** ************************************************************************
 *	Destructor for the Preconditioner class. 
 *  ************************************************************************ */
Preconditioner::~Preconditioner()
{
	ArrayUtils<double>::deltwotensor(vector);
}

/** ************************************************************************
 * The method to solve the system of equations associated with the
 * preconditioner.
 * 
 * Returns the value Solution class that is the solution to the
 * preconditioned system.  For the moment we just assume that the
 * preconditioner is the identity matrix. (This needs to be changed!)
 *
 * @param vector The Solution or right hand side of the system.
 * @return A Solution class member that is the solution to the
 *         preconditioned system.
 * ************************************************************************ */
Solution Preconditioner::solve(const Solution &current)
{
	Solution multiplied(current);
	int lupe;
	int innerlupe;
	for(lupe=0; lupe<getN(); lupe++)
	{
		for(innerlupe=0; innerlupe<Rank.COLNUMBER; innerlupe++)
		{
            multiplied(lupe,innerlupe) = current.getEntry(lupe,innerlupe);
		}
	}
	return(multiplied);
}


