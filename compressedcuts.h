#define MAXDIMENSIONS 5
#define MAX_ALLOWED_LEVELS 200

#define PTR_SIZE 4
#define HEADER_SIZE 4
#define BOUNDARY_SIZE 26
#define INTERVAL_SIZE 2

#define INTERNAL_NODE_SIZE (HEADER_SIZE + BOUNDARY_SIZE + PTR_SIZE)
#define LEAF_NODE_SIZE HEADER_SIZE

#define MAX_MEMBINS 4

#define TOO_MUCH 16

struct range{
  unsigned long long low;
  unsigned long long high;
};

struct pc_rule{
  int priority;
  struct range field[MAXDIMENSIONS];
  int siplen, diplen;
  unsigned sip[4], dip[4];
};

struct node 
{
  int depth;
  int problematic;
  int node_has_rule;
  pc_rule boundary;
  list <pc_rule*> classifier;
  list <node *> children;
  int cuts[MAXDIMENSIONS];
  // this is used only if this node is a result 
  // of cutting in 2D
  int Row;
  int Column;
  int Index;

  bool is_compressed;
};

struct TreeStat
{
  // independent vars
  int Id;
  int No_Rules;
  int Max_Depth;
  int Max_Levels;
  int Max_Cuts;
  int Rules_at_the_Leaf;
  int Rules_along_path;
  unsigned long long Total_Rule_Size;
  unsigned long long Total_Rules_Moved_Up;
  unsigned long long Total_Array_Size;
  unsigned long long Node_Count;
  unsigned long long Problematic_Node_Count;
  unsigned long long NonLeaf_Node_Count;
  unsigned long long Compressed_NonLeaf_Node_Count;
  unsigned long long Uncompressed_NonLeaf_Node_Count;
  map <unsigned,unsigned long long> interval_per_node;
  map <unsigned,unsigned long long> cuts_per_node;
  // dependent vars
  unsigned long long ruleptr_memory;
  unsigned long long array_memory;
  unsigned long long leaf_node_memory;
  unsigned long long compressed_int_node_memory;
  unsigned long long uncompressed_int_node_memory;
  unsigned long long total_memory;
  unsigned long long total_memory_in_KB;
};

struct MemBin
{
  int Max_Depth;
  int Max_Levels;
  unsigned long long total_memory;
  list<TreeStat*> Trees;
};

