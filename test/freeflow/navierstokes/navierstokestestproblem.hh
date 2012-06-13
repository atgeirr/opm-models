// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Christoph Grueninger                              *
 *   Copyright (C) 2009-2012 by Klaus Mosthaf                                *
 *   Copyright (C) 2010 by Katherina Baber                                   *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2012 by Andreas Lauser                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \brief  Definition of a simple Navier-Stokes problem
 */
#ifndef DUMUX_NAVIER_STOKES_TEST_PROBLEM_HH
#define DUMUX_NAVIER_STOKES_TEST_PROBLEM_HH

#include "n2constviscosity.hh"

#include <dumux/material/fluidsystems/h2on2fluidsystem.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/freeflow/stokes/stokesmodel.hh>

#if HAVE_ALUGRID
#include <dune/grid/alugrid/2d/alugrid.hh>
#elif HAVE_UG
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#else
#warning UG or ALUGrid necessary for this test.
#endif

#include <dune/grid/io/file/dgfparser.hh>

namespace Dumux {

template <class TypeTag>
class NavierStokesTestProblem;

// Specify the properties for the stokes problem
namespace Properties
{
NEW_TYPE_TAG(NavierStokesTestProblem, INHERITS_FROM(BoxStokes));

// Set the grid type
#if HAVE_ALUGRID
SET_TYPE_PROP(NavierStokesTestProblem, Grid, Dune::ALUCubeGrid<2,2>);
#elif HAVE_UG
SET_TYPE_PROP(NavierStokesTestProblem, Grid, Dune::UGGrid<2>);
#endif

// Set the problem property
SET_TYPE_PROP(NavierStokesTestProblem, Problem, Dumux::NavierStokesTestProblem<TypeTag>);

// Set calculation to Navier-Stokes, not Stokes
SET_BOOL_PROP(NavierStokesTestProblem, EnableNavierTerm, true);

SET_PROP(NavierStokesTestProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::GasPhase<Scalar, Dumux::N2ConstViscosity<Scalar> > type;
};

// Disable gravity
SET_BOOL_PROP(NavierStokesTestProblem, EnableGravity, false);

// Enable constraints
SET_BOOL_PROP(NavierStokesTestProblem, EnableConstraints, true);
}

/*!
 * \ingroup BoxStokesModel
 * \ingroup BoxTestProblems
 * \brief Stokes flow problem with modified nitrogen (N2) circulating in
 *        a cavity. (lid-driven cavity-flow)
 *
 * The example is taken from Ghia, Ghia, and Shin (1982), "High-Re solutions 
 * for incompressible flow using the Navier-Stokes equations and a multigrid
 * method", Journal of Computational Physics, Vol. 48, pp. 387-411.
 * 
 * The domain is two-dimensional and sized 1m times 1m. The boundary conditions 
 * for the momentum balances are Neumann zero boundary conditions except for
 * the top, which is floating from left to right with 1 m/s. The mass balance 
 * has outflow boundary conditions, which are replaced in the localresidual by 
 * the sum of the two momentum balances. All vertices at the bottom receive 
 * Dirichlet boundary conditions to set the pressure level.
 *
 * This problem uses the \ref BoxStokesModel with <code>EnableNavierStokes</code> 
 * set to <code>true</code>.
 */
template <class TypeTag>
class NavierStokesTestProblem : public StokesProblem<TypeTag>
{
    typedef StokesProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;

    enum {
        // Number of equations and grid dimension
        dimWorld = GridView::dimensionworld,

        // copy some indices for convenience
        pressureIdx = Indices::pressureIdx, 
        velocity0Idx = Indices::velocity0Idx,
        conti0EqIdx = Indices::conti0EqIdx,
        momentum0EqIdx = Indices::momentum0EqIdx
    };

    typedef typename GridView::ctype CoordScalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

public:
    NavierStokesTestProblem(TimeManager &timeManager)
        : ParentType(timeManager, GET_PROP_TYPE(TypeTag, GridCreator)::grid().leafView())
    { eps_ = 1e-6; }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    { return "navierstokes"; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a constant temperature of 10 degrees Celsius.
     */
    template <class Context>   
    Scalar temperature(const Context &context,
                       int spaceIdx, int timeIdx) const
    { return 273.15 + 10; } // -> 10 deg C

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Evaluate the boundary conditions.
     */
    template <class Context>
    void boundary(BoundaryRateVector &values, const Context &context, int spaceIdx, int timeIdx) const
    {
/*        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);      
        
        values.setOutflow(massBalanceIdx);
        values.setDirichlet(momentumXIdx);
        values.setDirichlet(momentumYIdx);
        // set pressure for all vertices at the bottom
        if (onLowerBoundary_(pos)) {
            values.setDirichlet(massBalanceIdx);
        }
*/
        values.setNoFlow(context, spaceIdx, timeIdx);
    }

    /*!
     * \brief Evaluate the constraints for a finite volume
     */
    template <class Context>
    void constraints(Constraints &values, const Context &context, int spaceIdx, int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);

        if (onUpperBoundary_(pos)) {
            // lid moves from left to right
            const Scalar lidVelocity = 1.0;
            values.setConstraint(momentum0EqIdx, velocity0Idx + 0, lidVelocity);
            values.setConstraint(momentum0EqIdx + 1, velocity0Idx + 1, 0);
            values.setConstraint(conti0EqIdx, pressureIdx, 1e5);
        }
    }
    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    template <class Context>
    void source(RateVector &values,
                const Context &context,
                int spaceIdx, int timeIdx) const
    { values = Scalar(0.0); }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    template <class Context>
    void initial(PrimaryVariables &priVars,
                 const Context &context,
                 int spaceIdx, int timeIdx) const
    { initial_(priVars); }

    // \}

private:
    // internal method for the initial condition
    void initial_(PrimaryVariables &priVars) const
    {
        priVars[pressureIdx] = 1e5;
        priVars[velocity0Idx + 0] = 0.0;
        priVars[velocity0Idx + 1] = 0.0;
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->bboxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->bboxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->bboxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->bboxMax()[1] - eps_;  }

    Scalar eps_;
};

} //end namespace

#endif