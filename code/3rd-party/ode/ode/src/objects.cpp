/*************************************************************************
 *                                                                       *
 * Open Dynamics Engine, Copyright (C) 2001,2002 Russell L. Smith.       *
 * All rights reserved.  Email: russ@q12.org   Web: www.q12.org          *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of EITHER:                                  *
 *   (1) The GNU Lesser General Public License as published by the Free  *
 *       Software Foundation; either version 2.1 of the License, or (at  *
 *       your option) any later version. The text of the GNU Lesser      *
 *       General Public License is included with this library in the     *
 *       file LICENSE.TXT.                                               *
 *   (2) The BSD-style license that is included with this library in     *
 *       the file LICENSE-BSD.TXT.                                       *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the files    *
 * LICENSE.TXT and LICENSE-BSD.TXT for more details.                     *
 *                                                                       *
 *************************************************************************/

// Object, body, and world methods.


#include <ode/common.h>
#include <ode/threading_impl.h>
#include <ode/objects.h>
#include <ode/config.h>
#include "matrix.h"
#include "objects.h"
#include "util.h"
#include "threading_impl.h"

static dThreadingImplementationID g_world_default_threading_impl = NULL;
static const dThreadingFunctionsInfo *g_world_default_threading_functions = NULL;

dxDampingParameters::dxDampingParameters(void *):
    linear_scale(REAL(0.0)),
    angular_scale(REAL(0.0)),
    linear_threshold(REAL(0.01) * REAL(0.01)),
    angular_threshold(REAL(0.01) * REAL(0.01))
{
}

bool dxWorld::InitializeDefaultThreading()
{
    dIASSERT(g_world_default_threading_impl == NULL);

    bool init_result = false;

    dThreadingImplementationID threading_impl = dThreadingAllocateSelfThreadedImplementation();

    if (threading_impl != NULL)
    {
        g_world_default_threading_functions = dThreadingImplementationGetFunctions(threading_impl);
        g_world_default_threading_impl = threading_impl;

        init_result = true;
    }

    return init_result;
}

void dxWorld::FinalizeDefaultThreading()
{
    dThreadingImplementationID threading_impl = g_world_default_threading_impl;

    if (threading_impl != NULL)
    {
        dThreadingFreeImplementation(threading_impl);

        g_world_default_threading_functions = NULL;
        g_world_default_threading_impl = NULL;
    }
}