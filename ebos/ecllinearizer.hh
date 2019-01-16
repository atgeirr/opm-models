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
 * \copydoc Ewoms::FvBaseLinearizer
 */
#ifndef EWOMS_ECL_LINEARIZER_HH
#define EWOMS_ECL_LINEARIZER_HH

#include <ewoms/disc/common/fvbaselinearizer.hh>

namespace Ewoms {

/*!
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief Add handling of focusTime to the linearizer.
 */
template<class TypeTag>
class EclLinearizer : public FvBaseLinearizer<TypeTag>
{
//! \cond SKIP_THIS
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef typename GET_PROP_TYPE(TypeTag, ThreadManager) ThreadManager;
    typedef typename GridView::template Codim<0>::Entity Element;

    typedef FvBaseLinearizer<TypeTag> Base;

public:
    EclLinearizer()
        : Base()
    {
    }


    void linearize(unsigned focusTimeIdx = 0)
    {
        linearizeDomain(focusTimeIdx);
        Base::linearizeAuxiliaryEquations();
    }

    void linearizeDomain(unsigned focusTimeIdx = 0)
    {
        focusTimeIndex_ = focusTimeIdx;
        Base::linearizeDomain();
    }


    // linearize an element in the interior of the process' grid partition
    void linearizeElement_(const Element& elem)
    {
        unsigned threadId = ThreadManager::threadId();

        ElementContext *elementCtx = this->elementCtx_[threadId];
        elementCtx->setFocusTimeIndex(this->focusTimeIndex_); // <---- the only line different from the base class
        auto& localLinearizer = this->model_().localLinearizer(threadId);

        // the actual work of linearization is done by the local linearizer class
        localLinearizer.linearize(*elementCtx, elem);

        // update the right hand side and the Jacobian matrix
        if (GET_PROP_VALUE(TypeTag, UseLinearizationLock))
            this->globalMatrixMutex_.lock();

        size_t numPrimaryDof = elementCtx->numPrimaryDof(/*timeIdx=*/0);
        for (unsigned primaryDofIdx = 0; primaryDofIdx < numPrimaryDof; ++ primaryDofIdx) {
            unsigned globI = elementCtx->globalSpaceIndex(/*spaceIdx=*/primaryDofIdx, /*timeIdx=*/0);

            // update the right hand side
            this->residual_[globI] += localLinearizer.residual(primaryDofIdx);

            // update the global Jacobian matrix
            for (unsigned dofIdx = 0; dofIdx < elementCtx->numDof(/*timeIdx=*/0); ++ dofIdx) {
                unsigned globJ = elementCtx->globalSpaceIndex(/*spaceIdx=*/dofIdx, /*timeIdx=*/0);

                this->jacobian_->addToBlock(globJ, globI, localLinearizer.jacobian(dofIdx, primaryDofIdx));
            }
        }

        if (GET_PROP_VALUE(TypeTag, UseLinearizationLock))
            this->globalMatrixMutex_.unlock();
    }

private:

    int focusTimeIndex_;
};

} // namespace Ewoms

#endif
