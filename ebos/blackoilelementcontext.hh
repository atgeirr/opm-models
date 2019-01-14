// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Ewoms::FvBaseElementContext
 */
#ifndef EWOMS_BLACKOIL_ELEMENT_CONTEXT_HH
#define EWOMS_BLACKOIL_ELEMENT_CONTEXT_HH

#include <ewoms/disc/common/fvbaseelementcontext.hh>

namespace Ewoms {

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief Adds the focusTimeIndex member to the base context.
 */
template<class TypeTag>
class BlackoilElementContext : public FvBaseElementContext<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef FvBaseElementContext<TypeTag> Base;
public:
    /*!
     * \brief The constructor.
     */
    explicit BlackoilElementContext(const Simulator& simulator)
        : Base(simulator)
    {
        focusTimeIdx_ = 0;
    }

    /*!
     * \brief Sets the time index on which the simulator is currently "focused" on
     *
     * I.e., in the case of automatic differentiation, all derivatives are with regard to
     * the primary variables of that time index. Only "primary" DOFs can be
     * focused on.
     */
    void setFocusTimeIndex(unsigned timeIdx)
    { focusTimeIdx_ = timeIdx; }

    /*!
     * \brief Returns the time index on which the simulator is currently "focused" on
     *
     * \copydetails setFocusDof()
     */ 
    unsigned focusTimeIndex() const
    { return focusTimeIdx_; }

public:
    void updateSingleIntQuants_(const PrimaryVariables& priVars, unsigned dofIdx, unsigned timeIdx)
    {
#ifndef NDEBUG
        if (Base::enableStorageCache_ && timeIdx != 0 && Base::problem().recycleFirstIterationStorage())
            throw std::logic_error("If caching of the storage term is enabled, only the intensive quantities "
                                   "for the most-recent substep (i.e. time index 0) are available!");
#endif
        Base::dofVars_[dofIdx].priVars[timeIdx] = priVars;
        Base::dofVars_[dofIdx].intensiveQuantities[timeIdx].update(*this, dofIdx, timeIdx, focusTimeIdx_);
    }

private:
    int focusTimeIdx_;	
};

} // namespace Ewoms

#endif
