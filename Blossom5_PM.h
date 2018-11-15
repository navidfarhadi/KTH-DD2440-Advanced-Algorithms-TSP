

#ifndef ASKHAKJSNTJAKSNBBAVASRA
#define ASKHAKJSNTJAKSNBBAVASRA

#ifdef _MSC_VER
#pragma warning(disable: 4311)
#pragma warning(disable: 4312)
#endif

#include <assert.h>
#include <string.h>
#include <limits.h>
#include "Blossom5_block.h"

#include "Blossom5_timer.h"
#include "Blossom5_PQ.h"
#include "Blossom5_LCA.h"


#define LCA_REPAIRS




#define IS_INT ( ((REAL)1 / 2) == 0 )
#define COST_FACTOR 2

#define PM_THRESHOLD ((REAL)1e-12)

template <typename FlowType, typename CostType> class MinCost
{
public:
	typedef int NodeId;
	typedef int EdgeId;

	MinCost(int NodeNum, int edgeNumMax, void (*err_function)(const char *) = NULL);

	
	~MinCost();

	void AddNodeExcess(NodeId i, FlowType excess);

	
	
	
	EdgeId AddEdge(NodeId i, NodeId j, FlowType cap, FlowType rev_cap, CostType cost);

	CostType Solve();

	

	FlowType GetRCap(EdgeId e);
	void SetRCap(EdgeId e, FlowType new_rcap);
	FlowType GetReverseRCap(EdgeId e);
	void SetReverseRCap(EdgeId e, FlowType new_rcap);
	void PushFlow(EdgeId e, FlowType delta);
	void UpdateCost(EdgeId e, FlowType cap_orig, CostType delta);

	CostType GetDual(NodeId i) { return nodes[i].pi; }

	
protected:
	

	struct Node;
	struct Arc;

	struct Node
	{
		Arc			*firstNonsaturated;
		Arc			*firstSaturated;

		Arc			*parent;
		Node		*next; 

		FlowType	excess;
		CostType	pi;
		int			flag;
		union
		{
			int		heap_ptr;
			Node*	next_permanent;
		};
#ifdef MINCOST_DEBUG
		int			id;
#endif
	};

	struct Arc
	{
		Node		*head;
		Arc			*prev;
		Arc			*next;
		Arc			*sister;	

		FlowType	r_cap;		
#ifdef MINCOST_DEBUG
		FlowType	cap_orig;
#endif
		CostType	cost;
		CostType GetRCost() { return cost + head->pi - sister->head->pi; }
	};

	int		nodeNum, edgeNum, edgeNumMax;
	Node	*nodes;
	Arc		*arcs;
	Node*	firstActive;
	int		counter;
	CostType cost;


	void	(*error_function)(const char *);										

	struct PriorityQueue
	{
		PriorityQueue();
		~PriorityQueue();
		void Reset();
		CostType GetKey(Node* i);
		void Add(Node* i, CostType key);
		void DecreaseKey(Node* i, CostType key);
		Node* RemoveMin(CostType& key);

	private:
		struct Item
		{
			Node*		i;
			CostType	key;
		}* array;
		int N, arraySize;
		void Swap(int k1, int k2);
	};

	PriorityQueue queue;

	

	void SetRCap(Arc* a, FlowType new_rcap);
	void PushFlow(Arc* a, FlowType delta);

	void Init();
	void DecreaseRCap(Arc* a, FlowType delta);
	void IncreaseRCap(Arc* a, FlowType delta);
	FlowType Augment(Node* start, Node* end);
	void Dijkstra(Node* start);

	void TestOptimality();
#ifdef MINCOST_DEBUG
	void TestCosts();
#endif
};


template <typename CostType> class DualMinCost : private MinCost<int, CostType>
{
public:
	typedef int NodeId;
	DualMinCost(int node_num, int constraint_num_max);
	~DualMinCost();

	void AddUnaryTerm(NodeId i, int objective_coef);
	void SetLowerBound(NodeId, CostType cmin);
	void SetUpperBound(NodeId, CostType cmax);
	void AddConstraint(NodeId i, NodeId j, CostType cmax); 

	void Solve();
	CostType GetSolution(NodeId i);

private:
	NodeId	source;
};

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::AddNodeExcess(NodeId _i, FlowType excess)
{
	assert(_i>=0 && _i<nodeNum);
	nodes[_i].excess += excess;
	if (nodes[_i].excess > 0 && !nodes[_i].next)
	{
		nodes[_i].next = firstActive;
		firstActive = &nodes[_i];
	}
}

template <typename FlowType, typename CostType> 
	inline typename MinCost<FlowType, CostType>::EdgeId MinCost<FlowType, CostType>::AddEdge(NodeId _i, NodeId _j, FlowType cap, FlowType rev_cap, CostType cost)
{
	assert(_i>=0 && _i<nodeNum);
	assert(_j>=0 && _j<nodeNum);
	assert(_i!=_j && edgeNum<edgeNumMax);
	assert(cap >= 0);
	assert(rev_cap >= 0);

	Arc *a = &arcs[2*edgeNum];
	Arc *a_rev = a+1;
	edgeNum ++;

	Node* i = nodes + _i;
	Node* j = nodes + _j;

	a -> sister = a_rev;
	a_rev -> sister = a;
	if (cap > 0)
	{
		if (i->firstNonsaturated) i->firstNonsaturated->prev = a;
		a -> next = i -> firstNonsaturated;
		i -> firstNonsaturated = a;
	}
	else
	{
		if (i->firstSaturated) i->firstSaturated->prev = a;
		a -> next = i -> firstSaturated;
		i -> firstSaturated = a;
	}
	a->prev = NULL;
	if (rev_cap > 0)
	{
		if (j->firstNonsaturated) j->firstNonsaturated->prev = a_rev;
		a_rev -> next = j -> firstNonsaturated;
		j -> firstNonsaturated = a_rev;
	}
	else
	{
		if (j->firstSaturated) j->firstSaturated->prev = a_rev;
		a_rev -> next = j -> firstSaturated;
		j -> firstSaturated = a_rev;
	}
	a_rev->prev = NULL;

	a -> head = j;
	a_rev -> head = i;
	a -> r_cap = cap;
	a_rev -> r_cap = rev_cap;
	a -> cost = cost;
	a_rev -> cost = -cost;
#ifdef MINCOST_DEBUG
	a->cap_orig = cap;
	a_rev->cap_orig = rev_cap;
#endif

	if (a->r_cap > 0 && a->GetRCost() < 0) PushFlow(a, a->r_cap);
	if (a_rev->r_cap > 0 && a_rev->GetRCost() < 0) PushFlow(a_rev, a_rev->r_cap);

	return edgeNum-1;
}


template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::DecreaseRCap(Arc* a, FlowType delta)
{
	a->r_cap -= delta;
	if (a->r_cap == 0)
	{
		Node* i = a->sister->head;
		if (a->next) a->next->prev = a->prev;
		if (a->prev) a->prev->next = a->next;
		else         i->firstNonsaturated = a->next;
		a->next = i->firstSaturated;
		if (a->next) a->next->prev = a;
		a->prev = NULL;
		i->firstSaturated = a;
	}
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::IncreaseRCap(Arc* a, FlowType delta)
{
	if (a->r_cap == 0)
	{
		Node* i = a->sister->head;
		if (a->next) a->next->prev = a->prev;
		if (a->prev) a->prev->next = a->next;
		else         i->firstSaturated = a->next;
		a->next = i->firstNonsaturated;
		if (a->next) a->next->prev = a;
		a->prev = NULL;
		i->firstNonsaturated = a;
	}
	a->r_cap += delta;
}

template <typename FlowType, typename CostType> 
	inline FlowType MinCost<FlowType, CostType>::GetRCap(EdgeId e)
{
	Arc* a = &arcs[2*e];
	return a->r_cap;
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::SetRCap(Arc* a, FlowType new_rcap)
{
	assert(new_rcap >= 0);
#ifdef MINCOST_DEBUG
	a->cap_orig += new_rcap - a->r_cap;
#endif
	if (a->r_cap == 0)
	{
		Node* i = a->sister->head;
		if (a->next) a->next->prev = a->prev;
		if (a->prev) a->prev->next = a->next;
		else         i->firstSaturated = a->next;
		a->next = i->firstNonsaturated;
		if (a->next) a->next->prev = a;
		a->prev = NULL;
		i->firstNonsaturated = a;
	}
	a->r_cap = new_rcap;
	if (a->r_cap == 0)
	{
		Node* i = a->sister->head;
		if (a->next) a->next->prev = a->prev;
		if (a->prev) a->prev->next = a->next;
		else         i->firstNonsaturated = a->next;
		a->next = i->firstSaturated;
		if (a->next) a->next->prev = a;
		a->prev = NULL;
		i->firstSaturated = a;
	}
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::SetRCap(EdgeId e, FlowType new_rcap)
{
	SetRCap(&arcs[2*e], new_rcap);
}

template <typename FlowType, typename CostType> 
	inline FlowType MinCost<FlowType, CostType>::GetReverseRCap(EdgeId e)
{
	Arc* a = &arcs[2*e+1];
	return a->r_cap;
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::SetReverseRCap(EdgeId e, FlowType new_rcap)
{
	SetRCap(&arcs[2*e+1], new_rcap);
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PushFlow(Arc* a, FlowType delta)
{
	if (delta < 0) { a = a->sister; delta = -delta; }
	DecreaseRCap(a, delta);
	IncreaseRCap(a->sister, delta);
	a->head->excess += delta;
	a->sister->head->excess -= delta;
	cost += delta*a->cost;
	if (a->head->excess > 0 && !a->head->next)
	{
		a->head->next = firstActive;
		firstActive = a->head;
	}
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PushFlow(EdgeId e, FlowType delta)
{
	PushFlow(&arcs[2*e], delta);
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::UpdateCost(EdgeId e, FlowType cap_orig, CostType delta)
{
	Arc* a = &arcs[2*e];
	cost += delta*(cap_orig-a->r_cap);
	a->cost += delta;
	a->sister->cost = -a->cost;

	if (a->GetRCost() > 0) a = a->sister;
	if (a->r_cap > 0 && a->GetRCost() < 0) PushFlow(a, a->r_cap);
}


template <typename FlowType, typename CostType> 
	inline MinCost<FlowType, CostType>::PriorityQueue::PriorityQueue()
{
	N = 0;
	arraySize = 16;
	array = (Item*) malloc(arraySize*sizeof(Item));
}

template <typename FlowType, typename CostType> 
	inline MinCost<FlowType, CostType>::PriorityQueue::~PriorityQueue()
{
	free(array);
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PriorityQueue::Reset()
{
	N = 0;
}

template <typename FlowType, typename CostType> 
	inline CostType MinCost<FlowType, CostType>::PriorityQueue::GetKey(Node* i)
{
	return array[i->heap_ptr].key;
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PriorityQueue::Swap(int k1, int k2)
{
	Item* a = array+k1;
	Item* b = array+k2;
	a->i->heap_ptr = k2;
	b->i->heap_ptr = k1;
	Node* i = a->i;   a->i   = b->i;   b->i   = i;
	CostType key = a->key; a->key = b->key; b->key = key;
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PriorityQueue::Add(Node* i, CostType key)
{
	if (N == arraySize)
	{
		arraySize *= 2;
		array = (Item*) realloc(array, arraySize*sizeof(Item));
	}
	int k = i->heap_ptr = N ++;
	array[k].i = i;
	array[k].key = key;
	while (k > 0)
	{
		int k2 = (k-1)/2;
		if (array[k2].key <= array[k].key) break;
		Swap(k, k2);
		k = k2;
	}
}

template <typename FlowType, typename CostType> 
	inline void MinCost<FlowType, CostType>::PriorityQueue::DecreaseKey(Node* i, CostType key)
{
	int k = i->heap_ptr;
	array[k].key = key;
	while (k > 0)
	{
		int k2 = (k-1)/2;
		if (array[k2].key <= array[k].key) break;
		Swap(k, k2);
		k = k2;
	}
}

template <typename FlowType, typename CostType> 
	inline typename MinCost<FlowType, CostType>::Node* MinCost<FlowType, CostType>::PriorityQueue::RemoveMin(CostType& key)
{
	if (N == 0) return NULL;

	Swap(0, N-1);
	N --;

	int k = 0;
	while ( 1 )
	{
		int k1 = 2*k + 1, k2 = k1 + 1;
		if (k1 >= N) break;
		int k_min = (k2 >= N || array[k1].key <= array[k2].key) ? k1 : k2;
		if (array[k].key <= array[k_min].key) break;
		Swap(k, k_min);
		k = k_min;
	}

	key = array[N].key;
	return array[N].i;
}

class PerfectMatching
{
public:

#ifdef PERFECT_MATCHING_DOUBLE
	typedef double REAL; 
	#define PM_INFTY ((REAL)1e100)
#else
	typedef int REAL;
	#define PM_INFTY (INT_MAX/2)
#endif

	typedef int NodeId;
	typedef int EdgeId;

	PerfectMatching(int nodeNum, int edgeNumMax);
	~PerfectMatching();

	
	EdgeId AddEdge(NodeId i, NodeId j, REAL cost);

	
	
	
	void Solve(bool finish=true);

	
	
	int GetSolution(EdgeId e); 
	NodeId GetMatch(NodeId i); 

	
	
	
	
	void GetDualSolution(int* blossom_parents, REAL* twice_y);
	int GetBlossomNum();

	
	
	
	
	void StartUpdate();
	void FinishUpdate();

	
	REAL GetTwiceSum(NodeId i); 
	EdgeId AddNewEdge(NodeId i, NodeId j, REAL cost, bool do_not_add_if_positive_slack=true); 
	void UpdateCost(EdgeId e, REAL delta_cost);


	
	
	
	
	struct Options
	{
		Options() : fractional_jumpstart(true),
		            dual_greedy_update_option(0),
		            dual_LP_threshold(0.00),
		            update_duals_before(false),
		            update_duals_after(false),
		            single_tree_threshold(1.00),
		            verbose(true)
		{}

		bool	fractional_jumpstart; 

		int 	dual_greedy_update_option; 
				                           
				                           

		double	dual_LP_threshold; 
		                           

		bool	update_duals_before; 
		bool	update_duals_after;  

		double	single_tree_threshold; 

		bool	verbose;
	} options;


	
	
	
	void Save(char* filename, int format=0); 

private:
	struct Node;
	struct Arc; 
	struct Edge; 
	struct Tree;
	struct TreeEdge;
	struct PQPointers;
	struct EdgeIterator;
	struct TreeEdgeIterator;
	struct LCATreeX;

	Node*	nodes;
	Edge*	edges;
	char*	edges_orig;
	DBlock<Node>* blossoms;
	Tree*	trees;
	DBlock<TreeEdge>* tree_edges;
	struct ExpandTmpItem
	{
		Node*	i;
		Node*	blossom_parent;
		Node*	blossom_grandparent;
	};
	Block<ExpandTmpItem>* expand_tmp_list; 

	int		node_num;
	int		edge_num, edge_num_max;
	int		tree_num, tree_num_max;

	Node*	removed_first;
	int		blossom_num;
	int		removed_num;

	void*	pq_buf;

	bool	first_solve;

	
	struct Stat
	{
		int		shrink_count;
		int		expand_count;
		int		grow_count;
		double	shrink_time;
		double	expand_time;
		double	dual_time;
	} stat;

	

	void InitGreedy(bool allocate_trees=true);

	void InitGlobal(); 
	Node* FindBlossomRootInit(Edge* a0);
	void ShrinkInit(Edge* a0, Node* tree_root);
	void ExpandInit(Node* b);
	void AugmentBranchInit(Node* i0, Node* tree_root);

	void Finish(); 

	void ProcessNegativeEdge(Edge* a);

	void GetRealEndpoints(Edge* a, Node*& tail, Node*& head);
	Node* FindBlossomRoot(Edge* a0);
	void Shrink(Edge* a0);
	void Expand(Node* b);
	void Augment(Edge* a0);
	void AugmentBranch(Node* i0);
	void GrowNode(Node* i);
	void GrowTree(Node* root, bool new_subtree);
	bool ProcessEdge00(Edge* a, bool update_boundary_edge=true); 
	void ProcessSelfloop(Node* b, Edge* a);

	void AddTreeEdge(Tree* t0, Tree* t1);

	void ComputeEpsSingle(); 
	void ComputeEpsCC(); 
	void ComputeEpsSCC(); 
	void ComputeEpsGlobal(); 
	bool UpdateDuals();

	void FreeRemoved();
	void CommitEps();

	void ReallocateEdges();

	void PrintAll();
};

int CheckPerfectMatchingOptimality(int node_num, int edge_num, int* edges, int* weights, PerfectMatching* pm, PerfectMatching::REAL threshold=(PerfectMatching::REAL)(1e-10));
double ComputePerfectMatchingCost(int node_num, int edge_num, int* edges, int* weights, PerfectMatching* pm);


struct PerfectMatching::Node
{
	unsigned int	is_outer : 1; 
	unsigned int	flag : 2; 
	unsigned int	is_tree_root : 1;
	unsigned int	is_processed : 1;
	unsigned int	is_blossom : 1;
	unsigned int	is_marked : 1;
	unsigned int	is_removed : 1;

	Edge*		first[2];
	union
	{
		Arc*	match; 
		Node*	blossom_grandparent;
	};
	REAL		y;

	union
	{
		struct 
		{
			Arc*	blossom_sibling;
			Node*	blossom_parent;
			union
			{
				Edge*	blossom_selfloops;
				Node*	blossom_ptr; 
#ifdef LCA_REPAIRS
				int		lca_preorder; 
#endif
			};
			REAL blossom_eps; 
			                  
		};
		struct 
		{
			union
			{
				struct 
				{
					Node*	first_tree_child;
					Node*	tree_sibling_prev; 
					Node*	tree_sibling_next;
				};
				Arc*	tree_parent; 
			};
			union
			{
				Tree*	tree;
				Edge*	best_edge;  
#ifdef LCA_REPAIRS
				int		lca_size; 
				LCATreeX*	lca;  
#endif
			};
		};
	};
};

struct PerfectMatching::Edge : PriorityQueue<REAL>::Item
{
	Node*	head[2];
	Node*	head0[2];
	Edge*	next[2];
	Edge*	prev[2];
};

typedef unsigned long POINTER_TYPE;
extern char dummy_array[2*(sizeof(void*)==sizeof(POINTER_TYPE))-1];

#define ARC_TO_EDGE_PTR(a)       ( (Edge*) ( ((POINTER_TYPE)(a)) & (~1)      ) )
#define ARC_TO_EDGE_DIR(a)       ( (int)   ( ((POINTER_TYPE)(a)) & 1         ) )
#define EDGE_DIR_TO_ARC(a, dir)  ( (Arc*)  ( (char*)(a) + (dir)) )

#define ARC_REV(a) ( (Arc*) ( ((POINTER_TYPE)(a)) ^ 1 ) )

#define ARC_TAIL(a)  (ARC_TO_EDGE_PTR(a)->head [1-ARC_TO_EDGE_DIR(a)])
#define ARC_TAIL0(a) (ARC_TO_EDGE_PTR(a)->head0[1-ARC_TO_EDGE_DIR(a)])
#define ARC_HEAD(a)  (ARC_TO_EDGE_PTR(a)->head [ARC_TO_EDGE_DIR(a)])
#define ARC_HEAD0(a) (ARC_TO_EDGE_PTR(a)->head0[ARC_TO_EDGE_DIR(a)])

struct PerfectMatching::PQPointers
{
	PriorityQueue<REAL> pq00; 
	union
	{
		PriorityQueue<REAL> pq01[2]; 
		struct 
		{
			PriorityQueue<REAL> pq0; 
			PriorityQueue<REAL> pq_blossoms;
		};
	};
};

struct PerfectMatching::Tree : PQPointers
{
	REAL		eps;
	TreeEdge*	first[2];
	Node*		root;

	PQPointers*	pq_current;
	int			dir_current;

	
	
	REAL		eps_delta;
	Tree*		next;
	union
	{
		int			id;
		TreeEdge*	dfs_parent;
	};
};

struct PerfectMatching::TreeEdge : PQPointers
{
	Tree*		head[2];
	TreeEdge*	next[2];
};





#define GET_PENULTIMATE_BLOSSOM(j)\
	{\
		Node* jtmp1 = j;\
		while ( 1 )\
		{\
			if (!j->blossom_grandparent->is_outer) j = j->blossom_grandparent;\
			else if (j->blossom_grandparent != j->blossom_parent) j->blossom_grandparent = j->blossom_parent;\
			else break;\
		}\
		Node* jtmp2;\
		for ( ; jtmp1!=j; jtmp1=jtmp2)\
		{\
			jtmp2 = jtmp1->blossom_grandparent;\
			jtmp1->blossom_grandparent = j;\
		}\
	}
#define GET_PENULTIMATE_BLOSSOM2(j)\
	{\
		Node* jtmp1 = j;\
		Node* jtmp_prev = NULL;\
		while ( 1 )\
		{\
			if (!j->blossom_grandparent->is_outer) { jtmp_prev = j; j = j->blossom_grandparent; }\
			else if (j->blossom_grandparent != j->blossom_parent) j->blossom_grandparent = j->blossom_parent;\
			else break;\
		}\
		if (jtmp_prev)\
		{\
			Node* jtmp2;\
			for ( ; jtmp1!=jtmp_prev; jtmp1=jtmp2)\
			{\
				jtmp2 = jtmp1->blossom_grandparent;\
				jtmp1->blossom_grandparent = jtmp_prev;\
			}\
		}\
	}




struct PerfectMatching::EdgeIterator
{
	Edge*	a_last;
	int		start_flag;
};

#define FOR_ALL_EDGES(i, a, dir, I)\
	for ( dir = (i->first[0]) ? 0 : 1, I.a_last = a = i->first[dir], I.start_flag = (a) ? 0 : 1;\
		a != I.a_last || (I.start_flag ++ == 0) || (dir ++ == 0 && (I.a_last = a = i->first[1]));\
		a = a->next[dir] )

#define CONTINUE_FOR_ALL_EDGES(i, a, dir, I)\
	for ( a = a->next[dir];\
		a != I.a_last || (I.start_flag ++ == 0) || (dir ++ == 0 && (I.a_last = a = i->first[1]));\
		a = a->next[dir] )

#define REMOVE_EDGE(i, a, dir)\
	{\
		if ((a)->prev[dir]==(a)) (i)->first[dir] = NULL;\
		else\
		{\
			(a)->prev[dir]->next[dir] = (a)->next[dir];\
			(a)->next[dir]->prev[dir] = (a)->prev[dir];\
			(i)->first[dir] = (a)->next[dir];\
		}\
	}

#define ADD_EDGE(i, a, dir)\
	{\
		if ((i)->first[dir])\
		{\
			(a)->prev[dir] = (i)->first[dir]->prev[dir];\
			(a)->next[dir] = (i)->first[dir];\
			(i)->first[dir]->prev[dir]->next[dir] = (a);\
			(i)->first[dir]->prev[dir] = (a);\
		}\
		else (i)->first[dir] = (a)->prev[dir] = (a)->next[dir] = (a);\
		(a)->head[1-(dir)] = (i);\
	}

#define MOVE_EDGE(i_old, i_new, a, dir)\
	{\
		REMOVE_EDGE(i_old, a, dir);\
		ADD_EDGE(i_new, a, dir);\
	}

#define GET_OUTER_HEAD(a, dir, j)\
	{\
		j = (a)->head[dir];\
		if (!j->is_outer)\
		{\
			Node* j_orig = j;\
			GET_PENULTIMATE_BLOSSOM(j);\
			j = j->blossom_parent;\
			int dir_rev = 1 - (dir);\
			MOVE_EDGE(j_orig, j, a, dir_rev);\
		}\
	}

#define GET_TREE_PARENT(child, parent)\
	{\
		Arc* a = (child)->tree_parent;\
		Edge* e = ARC_TO_EDGE_PTR(a);\
		int dir = ARC_TO_EDGE_DIR(a);\
		GET_OUTER_HEAD(e, dir, parent);\
	}


struct PerfectMatching::TreeEdgeIterator
{
	TreeEdge**	e_ptr;
};

#define FOR_ALL_TREE_EDGES(t, e, dir)\
	for ( dir = (t->first[0]) ? 0 : 1, e = t->first[dir];\
	      e || (dir ++ == 0 && (e = t->first[1]));\
	      e = e->next[dir] )

#define FOR_ALL_TREE_EDGES_X(t, e, dir, T)\
	for ( dir = (t->first[0]) ? 0 : 1, T.e_ptr = &t->first[dir], e = *T.e_ptr;\
	      e || (dir ++ == 0 && (e = *(T.e_ptr = &t->first[1])));\
	      e = *T.e_ptr )\
	if (e->head[dir] == NULL) { *T.e_ptr = e->next[dir]; tree_edges->Delete(e); }\
	else if ((T.e_ptr = &e->next[dir]))


#define MOVE_NODE_IN_TREE(i)\
	{\
		if ((i)->first_tree_child) (i) = (i)->first_tree_child;\
		else\
		{\
			while (!(i)->is_tree_root && !(i)->tree_sibling_next) { (i) = ARC_HEAD((i)->match); GET_TREE_PARENT(i, i); }\
			if ((i)->is_tree_root) break;\
			(i) = (i)->tree_sibling_next;\
		}\
	}

#define ADD_TREE_CHILD(i, j)\
	{\
		(j)->flag = 0;\
		(j)->tree = (i)->tree;\
		(j)->first_tree_child = NULL;\
		(j)->tree_sibling_next = (i)->first_tree_child;\
		if ((i)->first_tree_child)\
		{\
			(j)->tree_sibling_prev = (i)->first_tree_child->tree_sibling_prev;\
			(i)->first_tree_child->tree_sibling_prev = j;\
		}\
		else\
		{\
			(j)->tree_sibling_prev = j;\
		}\
		(i)->first_tree_child = j;\
	}


#define REMOVE_FROM_TREE(i)\
	{\
		if ((i)->tree_sibling_next) (i)->tree_sibling_next->tree_sibling_prev = (i)->tree_sibling_prev;\
		else\
		{\
			Node* i_NEXT = ARC_HEAD((i)->match); i_NEXT = ARC_HEAD(i_NEXT->tree_parent); i_NEXT = i_NEXT->first_tree_child;\
			i_NEXT->tree_sibling_prev = (i)->tree_sibling_prev;\
		}\
		if ((i)->tree_sibling_prev->tree_sibling_next) (i)->tree_sibling_prev->tree_sibling_next = (i)->tree_sibling_next;\
		else\
		{\
			Node* i_PARENT = ARC_HEAD((i)->match); i_PARENT = ARC_HEAD(i_PARENT->tree_parent);\
			i_PARENT->first_tree_child = (i)->tree_sibling_next;\
		}\
	}


#endif

