//////////////////////////////////////////////////////////////////////////////
// This file is part of the Maple Engine                              		//
//////////////////////////////////////////////////////////////////////////////
#include "Simplifier2.h"
#include <algorithm>
#include <cassert>
#include <functional>
#include <unordered_map>
#include <unordered_set>

namespace virtualGeometry
{
	class Heap
	{
		uint32_t  heapSize;
		uint32_t  numIndex;
		uint32_t *heap        = nullptr;
		float *   keys        = nullptr;
		uint32_t *heapIndexes = nullptr;

		inline auto pushUp(uint32_t i)
		{
			uint32_t idx = heap[i];
			uint32_t fa  = (i - 1) >> 1;
			while (i > 0 && keys[idx] < keys[heap[fa]])
			{
				heap[i]              = heap[fa];
				heapIndexes[heap[i]] = i;
				i = fa, fa = (i - 1) >> 1;
			}
			heap[i]              = idx;
			heapIndexes[heap[i]] = i;
		}

		inline auto pushDown(uint32_t i)
		{
			uint32_t idx = heap[i];
			uint32_t ls  = (i << 1) + 1;
			uint32_t rs  = ls + 1;
			while (ls < heapSize)
			{
				uint32_t t = ls;
				if (rs < heapSize && keys[heap[rs]] < keys[heap[ls]])
					t = rs;
				if (keys[heap[t]] < keys[idx])
				{
					heap[i]              = heap[t];
					heapIndexes[heap[i]] = i;
					i = t, ls = (i << 1) + 1, rs = ls + 1;
				}
				else
					break;
			}
			heap[i]              = idx;
			heapIndexes[heap[i]] = i;
		}

	  public:
		Heap()
		{
			heap = nullptr, keys = nullptr, heapIndexes = nullptr;
			heapSize = 0, numIndex = 0;
		}
		Heap(uint32_t InNumIndex)
		{
			heapSize    = 0;
			numIndex    = InNumIndex;
			heap        = new uint32_t[numIndex];
			keys        = new float[numIndex];
			heapIndexes = new uint32_t[numIndex];
			memset(heapIndexes, 0xff, numIndex * sizeof(uint32_t));
		}
		
		inline auto freeAll()
		{
			heapSize = 0, numIndex = 0;
			if (heap != nullptr)
				delete[] heap;
			if (keys != nullptr)
				delete[] keys;
			if (heapIndexes != nullptr)
				delete[] heapIndexes;
			heap = nullptr, keys = nullptr, heapIndexes = nullptr;
		}


		~Heap()
		{
			freeAll();
		}

		inline auto resize(uint32_t inNumIndex)
		{
			freeAll();
			heapSize    = 0;
			numIndex    = inNumIndex;
			heap        = new uint32_t[numIndex];
			keys        = new float[numIndex];
			heapIndexes = new uint32_t[numIndex];
			memset(heapIndexes, 0xff, numIndex * sizeof(uint32_t));
		}

		inline auto getKey(uint32_t idx) const
		{
			return keys[idx];
		}

		inline auto clear()
		{
			heapSize = 0;
			memset(heapIndexes, 0xff, numIndex * sizeof(uint32_t));
		}

		inline auto empty() const
		{
			return heapSize == 0;
		}

		inline auto isPresent(uint32_t idx) const
		{
			return heapIndexes[idx] != ~0u;
		}

		inline auto top() const
		{
			assert(heapSize > 0);
			return heap[0];
		}

		inline auto pop()
		{
			assert(heapSize > 0);

			uint32_t idx         = heap[0];
			heap[0]              = heap[--heapSize];
			heapIndexes[heap[0]] = 0;
			heapIndexes[idx]     = ~0u;
			pushDown(0);
		}

		inline auto add(float key, uint32_t idx)
		{
			assert(!isPresent(idx));

			uint32_t i       = heapSize++;
			heap[i]          = idx;
			keys[idx]        = key;
			heapIndexes[idx] = i;
			pushUp(i);
		}

		inline auto update(float key, uint32_t idx)
		{
			keys[idx]  = key;
			uint32_t i = heapIndexes[idx];
			if (i > 0 && key < keys[heap[(i - 1) >> 1]])
				pushUp(i);
			else
				pushDown(i);
		}

		inline auto remove(uint32_t idx)
		{
			assert(isPresent(idx));

			float    key = keys[idx];
			uint32_t i   = heapIndexes[idx];

			if (i == heapSize - 1)
			{
				--heapSize;
				heapIndexes[idx] = ~0u;
				return;
			}

			heap[i]              = heap[--heapSize];
			heapIndexes[heap[i]] = i;
			heapIndexes[idx]     = ~0u;
			if (key < keys[heap[i]])
				pushDown(i);
			else
				pushUp(i);
		}
	};

	struct vec3
	{
		float x, y, z;
	};

	struct dvec3
	{
		double x, y, z;

		dvec3()
		{
			x = 0, y = 0, z = 0;
		}

		dvec3(vec3 b)
		{
			x = b.x, y = b.y, z = b.z;
		}

		dvec3(double x, double y, double z) :
		    x(x),
		    y(y),
		    z(z)
		{
		}
	};

	template <typename T>
	inline void hashCode(std::size_t &seed, const T &v)
	{
		std::hash<T> hasher;
		seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
	}

	inline size_t getHash(const vec3 &v)
	{
		size_t seed = 0;
		hashCode(seed, v.x);
		hashCode(seed, v.y);
		hashCode(seed, v.z);
		return seed;
	}

	inline vec3 operator*(vec3 a, float b)
	{
		return vec3{a.x * b, a.y * b, a.z * b};
	}
	inline dvec3 operator*(dvec3 a, double b)
	{
		return dvec3{a.x * b, a.y * b, a.z * b};
	}

	inline vec3 operator+(vec3 a, vec3 b)
	{
		return vec3{a.x + b.x, a.y + b.y, a.z + b.z};
	}

	inline bool operator==(vec3 a, vec3 b)
	{
		return a.x == b.x && a.y == b.y && a.z == b.z;
	}

	inline bool operator==(dvec3 a, dvec3 b)
	{
		return a.x == b.x && a.y == b.y && a.z == b.z;
	}

	inline vec3 operator-(vec3 a, vec3 b)
	{
		return vec3{a.x - b.x, a.y - b.y, a.z - b.z};
	}
	inline dvec3 operator-(dvec3 a, dvec3 b)
	{
		return dvec3{a.x - b.x, a.y - b.y, a.z - b.z};
	}

	inline vec3 &operator+=(vec3 &a, vec3 b)
	{
		a.x += b.x, a.y += b.y, a.z += b.z;
		return a;
	}

	inline vec3 operator-(vec3 a)
	{
		return {-a.x, -a.y, -a.z};
	}

	template <typename T>
	struct vec4
	{
		T x, y, z, w;
	};

	template <typename T>
	struct mat4
	{
		union
		{
			T       m[16];
			vec4<T> cols[4];
		};
	};
}        // namespace virtualGeometry

namespace math
{
	using namespace virtualGeometry;

	inline vec3 cross(vec3 a, vec3 b)
	{
		return vec3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
	}

	inline dvec3 cross(dvec3 a, dvec3 b)
	{
		return dvec3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
	}

	inline float length(vec3 a)
	{
		return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
	}

	inline float length2(vec3 a)
	{
		return a.x * a.x + a.y * a.y + a.z * a.z;
	}

	inline vec3 normalize(vec3 a)
	{
		float rl = 1 / sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
		return a * rl;
	}

	inline dvec3 normalize(dvec3 a)
	{
		double rl = 1 / sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
		return a * rl;
	}

	inline float dot(vec3 a, vec3 b)
	{
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}

	inline double dot(dvec3 a, dvec3 b)
	{
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}

	template <typename T>
	inline bool invertColumnMajor(const T m[16], T invOut[16])
	{
		T       inv[16], det;
		int32_t i;

		inv[0]  = m[5] * m[10] * m[15] - m[5] * m[11] * m[14] - m[9] * m[6] * m[15] + m[9] * m[7] * m[14] + m[13] * m[6] * m[11] - m[13] * m[7] * m[10];
		inv[4]  = -m[4] * m[10] * m[15] + m[4] * m[11] * m[14] + m[8] * m[6] * m[15] - m[8] * m[7] * m[14] - m[12] * m[6] * m[11] + m[12] * m[7] * m[10];
		inv[8]  = m[4] * m[9] * m[15] - m[4] * m[11] * m[13] - m[8] * m[5] * m[15] + m[8] * m[7] * m[13] + m[12] * m[5] * m[11] - m[12] * m[7] * m[9];
		inv[12] = -m[4] * m[9] * m[14] + m[4] * m[10] * m[13] + m[8] * m[5] * m[14] - m[8] * m[6] * m[13] - m[12] * m[5] * m[10] + m[12] * m[6] * m[9];
		inv[1]  = -m[1] * m[10] * m[15] + m[1] * m[11] * m[14] + m[9] * m[2] * m[15] - m[9] * m[3] * m[14] - m[13] * m[2] * m[11] + m[13] * m[3] * m[10];
		inv[5]  = m[0] * m[10] * m[15] - m[0] * m[11] * m[14] - m[8] * m[2] * m[15] + m[8] * m[3] * m[14] + m[12] * m[2] * m[11] - m[12] * m[3] * m[10];
		inv[9]  = -m[0] * m[9] * m[15] + m[0] * m[11] * m[13] + m[8] * m[1] * m[15] - m[8] * m[3] * m[13] - m[12] * m[1] * m[11] + m[12] * m[3] * m[9];
		inv[13] = m[0] * m[9] * m[14] - m[0] * m[10] * m[13] - m[8] * m[1] * m[14] + m[8] * m[2] * m[13] + m[12] * m[1] * m[10] - m[12] * m[2] * m[9];
		inv[2]  = m[1] * m[6] * m[15] - m[1] * m[7] * m[14] - m[5] * m[2] * m[15] + m[5] * m[3] * m[14] + m[13] * m[2] * m[7] - m[13] * m[3] * m[6];
		inv[6]  = -m[0] * m[6] * m[15] + m[0] * m[7] * m[14] + m[4] * m[2] * m[15] - m[4] * m[3] * m[14] - m[12] * m[2] * m[7] + m[12] * m[3] * m[6];
		inv[10] = m[0] * m[5] * m[15] - m[0] * m[7] * m[13] - m[4] * m[1] * m[15] + m[4] * m[3] * m[13] + m[12] * m[1] * m[7] - m[12] * m[3] * m[5];
		inv[14] = -m[0] * m[5] * m[14] + m[0] * m[6] * m[13] + m[4] * m[1] * m[14] - m[4] * m[2] * m[13] - m[12] * m[1] * m[6] + m[12] * m[2] * m[5];
		inv[3]  = -m[1] * m[6] * m[11] + m[1] * m[7] * m[10] + m[5] * m[2] * m[11] - m[5] * m[3] * m[10] - m[9] * m[2] * m[7] + m[9] * m[3] * m[6];
		inv[7]  = m[0] * m[6] * m[11] - m[0] * m[7] * m[10] - m[4] * m[2] * m[11] + m[4] * m[3] * m[10] + m[8] * m[2] * m[7] - m[8] * m[3] * m[6];
		inv[11] = -m[0] * m[5] * m[11] + m[0] * m[7] * m[9] + m[4] * m[1] * m[11] - m[4] * m[3] * m[9] - m[8] * m[1] * m[7] + m[8] * m[3] * m[5];
		inv[15] = m[0] * m[5] * m[10] - m[0] * m[6] * m[9] - m[4] * m[1] * m[10] + m[4] * m[2] * m[9] + m[8] * m[1] * m[6] - m[8] * m[2] * m[5];

		det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

		if (det == 0)
			return false;

		det = 1.f / det;

		for (i = 0; i < 16; i++)
			invOut[i] = inv[i] * det;

		return true;
	}

	template <typename T>
	inline bool inverse(mat4<T> m, mat4<T> &inv)
	{
		return invertColumnMajor((T *) &m, (T *) &inv);
	}
}        // namespace math

namespace std
{
	template <>
	struct hash<virtualGeometry::vec3>
	{
		std::size_t operator()(const virtualGeometry::vec3 &vec) const
		{
			return virtualGeometry::getHash(vec);
		}
	};
};        // namespace std

namespace virtualGeometry
{
	struct Quadric
	{
		double a2, b2, c2, d2;
		double ab, ac, ad;
		double bc, bd, cd;

		Quadric()
		{
			memset(this, 0, sizeof(double) * 10);
		}

		Quadric(dvec3 p0, dvec3 p1, dvec3 p2)
		{
			dvec3  n = math::normalize(math::cross(p1 - p0, p2 - p0));
			double a = n.x;
			double b = n.y;
			double c = n.z;
			double d = -math::dot(n, p0);
			a2 = a * a, b2 = b * b, c2 = c * c, d2 = d * d;
			ab = a * b, ac = a * c, ad = a * d;
			bc = b * c, bd = b * d, cd = c * d;
		}

		inline auto add(const Quadric &b)
		{
			double *t1 = (double *) this;
			double *t2 = (double *) &b;
			for (uint32_t i = 0; i < 10; i++)
				t1[i] += t2[i];
		}

		inline auto get(vec3 &p)
		{
			mat4<double> m;
			mat4<double> m2;
			m.cols[0] = {a2, ab, ac, 0.0};
			m.cols[1] = {ab, b2, bc, 0};
			m.cols[2] = {ac, bc, c2, 0};
			m.cols[3] = {ad, bd, cd, 1};

			if (!math::inverse(m, m2))
				return false;

			vec4 v = m2.cols[3];

			p = {(float) v.x, (float) v.y, (float) v.z};
			return true;
		}

		inline auto evaluate(const vec3 &p)
		{
			double res = a2 * p.x * p.x + 2 * ab * p.x * p.y + 2 * ac * p.x * p.z + 2 * ad * p.x + b2 * p.y * p.y + 2 * bc * p.y * p.z + 2 * bd * p.y + c2 * p.z * p.z + 2 * cd * p.z + d2;
			return res <= 0.0 ? 0.0 : res;
		}
	};

	struct Simplifier2
	{
		Simplifier2(int8_t *data, uint32_t vertexCount, uint32_t *indices, uint32_t indexCount, uint32_t stride);
		auto lock(const vec3 &v) -> void;
		auto simplify(uint32_t target) -> float;
		auto addEdge(const vec3 &v0, const vec3 &v1, uint32_t index) -> bool;
		auto fixupTriangle(uint32_t triangleID) -> void;
		auto removeDuplicateVertex(uint32_t indexID) -> void;
		auto setVertexIndex(uint32_t indexID, uint32_t vertexId) -> void;
		auto isTriangleDuplicate(uint32_t triangleIdx) -> bool;
		auto compact() -> void;
		auto evaluate(const vec3 &v0, const vec3 &v1, bool merge) -> double;
		auto gatherAdjacentTriangles(const vec3 &v0, std::vector<uint32_t> &tris, bool &lock) -> void;
		auto merge(const std::vector<vec3> &vertex, const std::function<void(const std::vector<uint32_t> &)> &callback) -> void;
		auto getVertex(uint32_t i) -> vec3 &;

		int8_t *                                               data        = nullptr;
		uint32_t                                               vertexCount = 0;
		uint32_t *                                             indices     = nullptr;
		uint32_t                                               indexCount  = 0;
		std::vector<bool>                                      deletedMarks;        //triangle marks..
		std::vector<uint8_t>                                   flags;
		std::vector<Quadric>                                   errorQ;
		std::vector<std::pair<vec3, vec3>>                     edges;
		std::vector<uint32_t>                                  reevaluateEdge;
		std::unordered_map<vec3, std::unordered_set<uint32_t>> edgeToID0;
		std::unordered_map<vec3, std::unordered_set<uint32_t>> edgeToID1;
		std::vector<uint32_t>                                  vertexRefCount;
		std::unordered_map<vec3, std::unordered_set<int32_t>>  vertexToID;        //there are duplicated vertex....
		std::unordered_map<vec3, std::unordered_set<int32_t>>  indexIDs;          // Vertex mapping to indexBuffer's index
		uint32_t                                               remainingVertex;
		uint32_t                                               remainingTriangles;
		uint32_t                                               triangles = 0;
		float                                                  maxError  = 0.f;
		Heap                                                   edgeHeap;
		uint32_t                                               stride = 0;
	};

	namespace
	{
#define BIT(x) (1 << x)
		enum Flags
		{
			NONE     = BIT(0),
			ADJACENT = BIT(1),
			LOCK     = BIT(2)
		};

		inline auto getNextEdge(uint32_t i) -> uint32_t
		{
			uint32_t mod3 = i % 3;
			return i - mod3 + ((1 << mod3) & 3);
		}
		inline auto getNextEdge(uint32_t i, uint32_t offset) -> uint32_t
		{
			return i - i % 3 + (i + offset) % 3;
		}

		inline auto isValid(const vec3 &p0, const vec3 &p1, const vec3 &p)
		{
			return math::length(p - p0) + math::length(p - p1) <= 2 * math::length(p0 - p1);
		}
	}        // namespace

	Simplifier2::Simplifier2(int8_t *data, uint32_t vertexCount, uint32_t *indices, uint32_t indexCount, uint32_t stride) :
	    stride(stride), data(data), vertexCount(vertexCount), indexCount(indexCount), indices(indices)
	{
		triangles = indexCount / 3;
		flags.resize(indexCount, 0);
		deletedMarks.resize(triangles, false);
		vertexRefCount.resize(vertexCount, 0);

		remainingVertex    = vertexCount;
		remainingTriangles = triangles;

		for (uint32_t i = 0; i < vertexCount; i++)
		{
			vertexToID[getVertex(i)].emplace(i);
		}

		//edgeConstruct
		for (auto idx = 0; idx < indexCount; idx++)
		{
			auto vertexID = indices[idx];
			vertexRefCount[vertexID]++;
			const auto &v = getVertex(vertexID);

			indexIDs[v].emplace(idx);

			auto index2 = getNextEdge(idx);
			auto p0     = v;
			auto p1     = getVertex(indices[index2]);

			if (getHash(p0) > getHash(p1))
				std::swap(p0, p1);

			if (addEdge(p0, p1, edges.size()))
				edges.push_back({p0, p1});
		}
	}

	auto Simplifier2::getVertex(uint32_t i) -> vec3 &
	{
		vec3 *newData = reinterpret_cast<vec3 *>(data + i * stride);
		return newData[0];
	}

	auto Simplifier2::addEdge(const vec3 &v0, const vec3 &v1, uint32_t index) -> bool
	{
		for (uint32_t i : edgeToID0[v0])
		{
			auto &e = edges[i];
			if (e.first == v0 && e.second == v1)
				return false;
		}
		edgeToID0[v0].emplace(index);
		edgeToID1[v1].emplace(index);
		return true;
	}

	auto Simplifier2::lock(const vec3 &v) -> void
	{
		for (auto &id : indexIDs[v])
		{
			if (getVertex(indices[id]) == v)
			{
				flags[id] = Flags::LOCK;
			}
		}
	}

	auto Simplifier2::simplify(uint32_t target) -> float
	{
		errorQ.resize(triangles);

		for (uint32_t i = 0; i < triangles; i++)
			fixupTriangle(i);

		if (remainingTriangles <= target)
		{
			compact();
			/*vertices.resize(remainingVertex);
				indices.resize(remainingTriangles * 3);*/
			return maxError;
		}

		edgeHeap.resize(edges.size());
		uint32_t i = 0;
		for (auto &e : edges)
		{
			double error = evaluate(e.first, e.second, false);
			edgeHeap.add(error, i);
			i++;
		}

		maxError = 0.f;
		while (!edgeHeap.empty())
		{
			auto  index = edgeHeap.top();
			float key   = edgeHeap.getKey(index);
			if (key >= 1e6)
				break;

			edgeHeap.pop();

			auto &e = edges[index];

			edgeToID0[e.first].erase(index);
			edgeToID1[e.second].erase(index);

			double error = evaluate(e.first, e.second, true);
			if (error > maxError)
				maxError = error;

			if (remainingTriangles <= target)
				break;

			for (auto i : reevaluateEdge)
			{
				auto & e     = edges[i];
				double error = evaluate(e.first, e.second, false);
				edgeHeap.add(error, i);
			}
			reevaluateEdge.clear();
		}
		compact();

		/*	vertices.resize(remainingVertex);
			indices.resize(remainingTriangles * 3);*/

		return maxError;
	}

	auto Simplifier2::fixupTriangle(uint32_t triangleID) -> void
	{
		assert(!deletedMarks[triangleID], "triangle was already removed...");

		const auto &p0 = getVertex(indices[triangleID * 3]);
		const auto &p1 = getVertex(indices[triangleID * 3 + 1]);
		const auto &p2 = getVertex(indices[triangleID * 3 + 2]);

		bool removed = false;
		if (!removed)
		{
			removed = (p0 == p1) || (p1 == p2) || (p2 == p0);
		}

		if (!removed)
		{
			for (auto k = 0; k < 3; k++)
			{
				removeDuplicateVertex(triangleID * 3 + k);
			}
			removed = isTriangleDuplicate(triangleID);
		}

		if (removed)
		{
			deletedMarks[triangleID] = true;
			remainingTriangles--;
			for (uint32_t k = 0; k < 3; k++)
			{
				auto corner = triangleID * 3 + k;
				auto index  = indices[corner];
				indexIDs[getVertex(index)].erase(corner);
				setVertexIndex(corner, ~0u);
			}
		}
		else
		{
			errorQ[triangleID] = Quadric(p0, p1, p2);
		}
	}

	auto Simplifier2::removeDuplicateVertex(uint32_t indexID) -> void
	{
		const auto  index = indices[indexID];
		const auto &v     = getVertex(index);

		for (auto vertexId : vertexToID[v])
		{
			if (vertexId == index)
				break;

			if (v == getVertex(vertexId))
			{
				setVertexIndex(indexID, vertexId);
				break;
			}
		}
	}

	auto Simplifier2::setVertexIndex(uint32_t indexID, uint32_t vertexId) -> void
	{
		auto &verIdx = indices[indexID];
		assert(verIdx != ~0u, "Index Error");
		assert(vertexRefCount[verIdx] > 0, "have not been referenced..");

		if (verIdx == vertexId)
			return;
		if (--vertexRefCount[verIdx] == 0)
		{
			vertexToID[getVertex(verIdx)].erase(verIdx);
			remainingVertex--;
		}

		verIdx = vertexId;
		if (verIdx != ~0u)
			vertexRefCount[verIdx]++;
	}

	auto Simplifier2::isTriangleDuplicate(uint32_t triangleIdx) -> bool
	{
		uint32_t i0 = indices[triangleIdx * 3 + 0];
		uint32_t i1 = indices[triangleIdx * 3 + 1];
		uint32_t i2 = indices[triangleIdx * 3 + 2];

		for (auto i : indexIDs[getVertex(i0)])
		{
			if (i != triangleIdx * 3)
			{
				if (i0 == indices[i] && i1 == indices[getNextEdge(i)] && i2 == indices[(getNextEdge(i, 2))])
					return true;
			}
		}
		return false;
	}

	auto Simplifier2::compact() -> void
	{
		uint32_t count = 0;
		for (uint32_t i = 0; i < vertexCount; i++)
		{
			if (vertexRefCount[i] > 0)
			{
				if (i != count)
					getVertex(count) = getVertex(i);
				vertexRefCount[i] = count++;
			}
		}

		assert(count == remainingVertex, "vertex count did not match..");

		uint32_t triangleCount = 0;
		for (uint32_t i = 0; i < triangles; i++)
		{
			if (!deletedMarks[i])
			{
				for (uint32_t k = 0; k < 3; k++)
				{
					indices[triangleCount * 3 + k] = vertexRefCount[indices[i * 3 + k]];
				}
				triangleCount++;
			}
		}
		assert(triangleCount == remainingTriangles, "triangle count did not match ...");
	}

	auto Simplifier2::evaluate(const vec3 &v0, const vec3 &v1, bool needMerge) -> double
	{
		if (v0 == v1)
			return 0.f;

		double error = 0;

		std::vector<uint32_t> adjacentTriangles;

		bool lock0 = false;
		bool lock1 = false;

		gatherAdjacentTriangles(v0, adjacentTriangles, lock0);
		gatherAdjacentTriangles(v1, adjacentTriangles, lock1);

		if (adjacentTriangles.empty())
			return 0.f;

		/*	if (adjacentTriangles.size() > 24)
			{
				error += 0.5 * (adjacentTriangles.size() - 24);
			}*/

		/*	glm::dmat4 quadric(0);

				for (auto idx : adjacentTriangles)
				{
					quadric += errorQ[idx];
				}*/

		Quadric q;
		for (uint32_t i : adjacentTriangles)
		{
			q.add(errorQ[i]);
		}

		vec3 p = (v0 + v1) * 0.5f;
		if (lock0 && lock1)
			error += 1e8;
		if (lock0 && !lock1)
			p = v0;
		else if (!lock0 && lock1)
			p = v1;
		else if (!q.get(p))
			p = (v0 + v1) * 0.5f;

		if (!isValid(v0, v1, p))
		{
			p = (v0 + v1) * 0.5f;
		}

		error += q.evaluate(p);
		//quadricEvaluate(p, quadric);

		if (needMerge)
		{
			merge({v0, v1}, [&](const std::vector<uint32_t> &moveEdgeQ) {
				for (auto i : adjacentTriangles)
				{
					for (auto k = 0; k < 3; k++)
					{
						uint32_t index = i * 3 + k;
						auto &   pos   = getVertex(indices[index]);
						if (pos == v0 || pos == v1)
						{
							pos = p;
							if (lock0 || lock1)
								flags[index] |= LOCK;
						}
					}
				}

				for (auto i : moveEdgeQ)
				{
					auto &e = edges[i];

					if (e.first == v0 || e.first == v1)
						e.first = p;

					if (e.second == v0 || e.second == v1)
						e.second = p;
				}
			});

			std::vector<uint32_t> adjacentVertices;
			for (uint32_t i : adjacentTriangles)
			{
				for (uint32_t k = 0; k < 3; k++)
				{
					adjacentVertices.push_back(indices[i * 3 + k]);
				}
			}

			std::sort(adjacentVertices.begin(), adjacentVertices.end());
			adjacentVertices.erase(std::unique(adjacentVertices.begin(), adjacentVertices.end()), adjacentVertices.end());

			for (uint32_t vertexId : adjacentVertices)
			{
				for (uint32_t i : edgeToID0[getVertex(vertexId)])
				{
					if (edges[i].first == getVertex(vertexId))
					{
						if (edgeHeap.isPresent(i))
						{
							edgeHeap.remove(i);
							reevaluateEdge.push_back(i);
						}
					}
				}
				for (uint32_t i : edgeToID1[getVertex(vertexId)])
				{
					if (edges[i].second == getVertex(vertexId))
					{
						if (edgeHeap.isPresent(i))
						{
							edgeHeap.remove(i);
							reevaluateEdge.push_back(i);
						}
					}
				}
			}

			for (uint32_t i : adjacentTriangles)
			{
				fixupTriangle(i);
			}
		}

		for (auto i : adjacentTriangles)
		{
			flags[i * 3] &= (~ADJACENT);
		}
		return error;
	}

	auto Simplifier2::gatherAdjacentTriangles(const vec3 &v0, std::vector<uint32_t> &tris, bool &lock) -> void
	{
		for (uint32_t i : indexIDs[v0])
		{
			if (getVertex(indices[i]) == v0)
			{
				uint32_t triangleIdx = i / 3;
				if ((flags[triangleIdx * 3] & ADJACENT) == 0)
				{
					flags[triangleIdx * 3] |= ADJACENT;
					tris.push_back(triangleIdx);
				}

				lock = (flags[i] & LOCK) == LOCK;
			}
		}
	}

	auto Simplifier2::merge(const std::vector<vec3> &vertex, const std::function<void(const std::vector<uint32_t> &)> &callback) -> void
	{
		std::vector<uint32_t> moveVertexQ;
		std::vector<uint32_t> moveIndexQ;
		std::vector<uint32_t> moveEdgeQ;

		//begin merge...

		for (auto &v : vertex)
		{
			{
				auto ids = vertexToID[v];
				for (auto id : ids)
				{
					if (getVertex(id) == v)
					{
						vertexToID[v].erase(id);
						moveVertexQ.emplace_back(id);
					}
				}
			}

			{
				auto ids = indexIDs[v];
				for (auto id : ids)
				{
					if (getVertex(indices[id]) == v)
					{
						indexIDs[v].erase(id);
						moveIndexQ.emplace_back(id);
					}
				}
			}

			{
				auto ids = edgeToID0[v];
				for (auto id : ids)
				{
					if (edges[id].first == v)
					{
						edgeToID0[edges[id].first].erase(id);
						edgeToID1[edges[id].second].erase(id);
						moveEdgeQ.push_back(id);
					}
				}
			}

			{
				auto ids = edgeToID1[v];
				for (auto id : ids)
				{
					if (edges[id].second == v)
					{
						edgeToID0[edges[id].first].erase(id);
						edgeToID1[edges[id].second].erase(id);
						moveEdgeQ.push_back(id);
					}
				}
			}
		}

		callback(moveEdgeQ);

		//end
		for (auto i : moveVertexQ)
		{
			vertexToID[getVertex(i)].emplace(i);
		}

		for (auto i : moveIndexQ)
		{
			indexIDs[getVertex(indices[i])].emplace(i);
		}

		for (auto i : moveEdgeQ)
		{
			auto &e = edges[i];
			if (e.first == e.second || !addEdge(e.first, e.second, i))
			{
				edgeHeap.remove(i);
			}
		}
	}

}        // namespace virtualGeometry

namespace virtualGeometry
{
	float simplifyByQem(uint32_t target, void *data, uint32_t vertexCount, uint32_t *indices, uint32_t indexCount, uint32_t stride,
	                    const float *lockedVertices, uint32_t count,
	                    uint32_t &outIndex, uint32_t &outVertex)
	{
		Simplifier2 simplifier(reinterpret_cast<int8_t *>(data), vertexCount, indices, indexCount, stride);

		auto locked = reinterpret_cast<const vec3 *>(lockedVertices);

		for (auto i = 0; i < count; i++)
			simplifier.lock(locked[i]);

		float error = simplifier.simplify(target);
		outIndex    = simplifier.remainingTriangles * 3;
		outVertex   = simplifier.remainingVertex;
		return std::sqrt(error);
	}        // namespace virtualGeometry
}        // namespace virtualGeometry
