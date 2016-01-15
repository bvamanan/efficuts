#!/usr/local/bin/perl
use bignum;
use POSIX qw(ceil floor);

$logdir = $ARGV[0];

# bunch of constants
$FIELDS_IN_TREE           = 10;

$BOUNDARY_SIZE            = 16;
$HEADER_SIZE              = 2;
$LEAF_HEADER_SIZE         = 1;
$PTR_SIZE                 = 4;
$RULE_ENTRY_SIZE          = 13;

$ARRAY_ENTRY_SIZE         = $BOUNDARY_SIZE + $HEADER_SIZE + $PTR_SIZE;

$RULES_PER_ARRAY_ENTRY    = ($ARRAY_ENTRY_SIZE - $PTR_SIZE - $LEAF_HEADER_SIZE)/ $RULE_ENTRY_SIZE;

$BUS_WIDTH                = $ARRAY_ENTRY_SIZE;

$output = `ls $logdir/*.out`;
@files = split(/\s+/,$output);

for $file (@files)
{
  # Infer stuff from file name

  @list_1 = split(/\//,$file);

  @list_2 = split(/_/,$list_1[$#list_1]);

  $category = $list_2[0];

  @list_3 = split(/\./,$list_2[2]);

  $numrules = int($list_3[0]);
  
  $hypercuts = int($list_2[3]);

  $compressionON = int($list_2[4]);

  @list_4 = split(/\./,$list_2[5]);

  $num_intervals = int($list_4[0]);

  #reset all data read previosly
  reset 'Tree';
  reset 'balajee_';
  reset 'memory_';

  open INFILE, "$file" or die $!;
  while (<INFILE>)
  {
    if (/^hypercuts = (\d+)$/) { 
      if ($hypercuts != int($1)) { die "$file: hypercuts mismatch in $file($hypercuts != $1)\n"; }
    }

    if (/^compressionON = (\d+)$/) {
      if ($compressionON != int($1)) { die "$file: compressionON mismatch in $file\n"; }
    }

    if (/^num_intervals = (\d+)$/) {
      $balajee_num_intervals = int($1);
      if ($num_intervals != int($1)) { die "$file: num_intervals mismatch in $file\n"; }
    }

    if (/^number of rules read from file = (\d+)$/) 
    {
      $balajee_in_rules = int($1);
      if ($1 < (0.40 * $numrules)) { warn "$file: rules read = $1!\n"; }
    }

    if (/^Tree: (\d+)$/) 
    { 
      reset 'tree_';
      $tree_label = $1; 
    }

    if (/^Rules: (\d+)$/) { push(@tree_record,$1); }

    if (/^Levels: (\d+)$/) { push(@tree_record,$1); }

    if (/^Rules_at_the_Leaf: (\d+)$/) { push(@tree_record,$1); }

    if (/^Rules_along_path: (\d+)$/) 
    { 
      push(@tree_record,$1); 
      if (int($1) != 0 && $hypercuts == 1 && $compressionON == 1)
      {
        die "hypercuts = $hypercuts;compressionON = $compressionON;Rules_along_path = $1\n";
      }
    }

    if (/^Total_Rule_Size: (\d+)$/) { push(@tree_record,$1); }

    if (/^Total_Array_Size: (\d+)$/) { push(@tree_record,$1); }

    if (/^Node_Count: (\d+)$/) { push(@tree_record,$1); }

    if (/^NonLeaf_Node_Count: (\d+)$/) { push(@tree_record,$1); }

    if (/^Compressed_NonLeaf_Node_Count: (\d+)$/) { push(@tree_record,$1); }

    if (/^Uncompressed_NonLeaf_Node_Count: (\d+)$/) { 
      
      push(@tree_record,$1); 

      if (@tree_record != $FIELDS_IN_TREE)
      {
        die "size of record != $FIELDS_IN_TREE(@tree_record)!\n";
      }

      if (exists $TreeData{$tree_label})
      {
        die "Multiple trees with the same tag!\n";
      }

      $TreeData{$tree_label} = [ @tree_record ]; 

    }

  }
  close INFILE;

  # print ruleset,size
  printf("%5s,%7s,",$category,$numrules);

  # Aggregate information from all trees
  for $key (sort {$a <=> $b} keys %TreeData)
  {
    # count the rules in each tree 
    $balajee_num_rules += $TreeData{$key}[0];

    # memory accesses = (levels * array_entry_size + (rules_at_leaf - rules_per_array_entry) * rule_size)/bus_width

    $memory_accesses_per_tree = ceil( ( $TreeData{$key}[1] * $ARRAY_ENTRY_SIZE + 
                                      ( $TreeData{$key}[2] + $TreeData{$key}[3] - $RULES_PER_ARRAY_ENTRY )
                                        * $RULE_ENTRY_SIZE ) / $BUS_WIDTH );

    if ($memory_accesses_per_tree <= 0) { die "$file: [$key] memory accesses <= 0\n"; }

    $memory_accesses += $memory_accesses_per_tree;

    # memory size = (array_entries * array_entry_size) + 
    #               (rules - (leaves * rules_per_array_entry)) * rule_size
    
    # $memory_arrays_per_tree = $TreeData{$key}[5] * $ARRAY_ENTRY_SIZE;
    $memory_arrays_per_tree = ($TreeData{$key}[5] - $TreeData{$key}[8] * ($balajee_num_intervals + 1)) * $PTR_SIZE
                            + ($TreeData{$key}[8] * ($balajee_num_intervals + 1)) * $ARRAY_ENTRY_SIZE 
                            + ($TreeData{$key}[9] + ($TreeData{$key}[6] - $TreeData{$key}[7])) * $ARRAY_ENTRY_SIZE;

    $memory_rules_per_tree  = ($TreeData{$key}[4] - ($TreeData{$key}[6] - $TreeData{$key}[7]) * $RULES_PER_ARRAY_ENTRY)
                                  * $RULE_ENTRY_SIZE;

    if ($memory_rules_per_tree < 0) { $memory_rules_per_tree = 0; }

    $memory_per_tree = $memory_arrays_per_tree + $memory_rules_per_tree;

    if ($memory_per_tree < $TreeData{$key}[0] * $RULE_ENTRY_SIZE)
    {
      warn "$file Tree[$key]: memory of tree($memory_per_tree) < $TreeData{$key}[0] * $RULE_ENTRY_SIZE!\n";
    }

    $memory_arrays += $memory_arrays_per_tree;

    $memory_rules += $memory_rules_per_tree;

  }

  if ($balajee_num_rules != $balajee_in_rules)
  {
    die "$file: Incoming number of rules = sum(rules_in_tree); $balajee_num_rules != $balajee_in_rules\n";
  }

  $memory_size = ($memory_arrays + $memory_rules)/(1024 * 1024);
  
  printf("%8d,%12.2f,%5d,",$balajee_in_rules,$memory_size,$memory_accesses);

  printf("\n");

} 



