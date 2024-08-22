//////////////////////////////////////////////////////////////////////////////
// This file is part of the Maple Engine                              		//
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include <cstdint>

namespace virtualGeometry
{
	float simplifyByQem(uint32_t target, void *data, uint32_t vertexCount, uint32_t *indices, uint32_t indexCount, uint32_t stride, 
		const float *lockedVertices, uint32_t count,
		uint32_t &outIndex, uint32_t &outVertex);
}        // namespace virtual_geometry
