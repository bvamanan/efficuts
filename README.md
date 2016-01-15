# efficuts
http://dl.acm.org/citation.cfm?id=1851208
The code & post-processing script are attached. 

After building, the run syntax is:
./compressedcuts -b16 -s8 -f7 -r acl1_seed_1000.filter -m1 -c1 -n0.5 -i0.05 -t0 -u0 -g1 -z1

The explanation of each of the arguments is in the source file. It covers HiCuts, HyperCuts and EffiCuts based on arguments. 

If I recall correctly, 

m - controls if you want to cut in more than one dimension 
    m0 - HiCuts
    m1 - HyperCuts/EffiCuts

g,z,c   - control EffiCuts optimizations
          c0,g0,z0 - turns off ALL EffiCuts' techniques
          c1,g1,z1 - EffiCuts 

So, m0,c0,g0,z0 = HiCuts
    m1,c0,g0,z0 = HyperCuts
    m1,c1,g1,z1 = EffiCuts

These are immediately apparent from the source code. 

We used ClassBench generated rulesets for our studies. You can download and use the same.
