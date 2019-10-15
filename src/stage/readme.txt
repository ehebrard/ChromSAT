== MAN BNP =====================================================================

bnp <filename> [-i  <input format>   || -o  <output format>  ||
                -sm <selector mode>  || -gm <generator mode> ||
                -d                   || -mono]

with :
<input format>   among tgf, dot, clq and dimacs
<output format>  among dot, sol and std
<selector mode>  among swap, sort and near
<generator mode> among cplex and greedy
-d               limits the output printed on the shell
-mono            force cplex to use only 1 thread

== CREATE A SOLVER =============================================================

>> How to create a solver?
Simply create an object without any call to a constructor.
Most of its attributes will be created when a graph is loaded.

== LOAD A GRAPH ================================================================

>> How to load a graph?
Call BnP::load(String filename, InFormat in) with:
#-------------#----------------------------------------------------------------#
| filename    | Simply the name of the file                                    |
#-------------#----------------------------------------------------------------#
| in          | I_TGF, I_DOT, I_CLQ or I_DIMACS                                |
#-------------#----------------------------------------------------------------#
/!\ in contains all the planed format. But every single reader listed are    /!\
/!\ not coded yet. The main format DIMACS do work.                           /!\

== TUNE THE SOLVER =============================================================

>> How to use your own column-generating methods?
Your methods must have the Generator form: (IloNumArray, ms_graph)-> vector<int>
with:
#-------------#----------------------------------------------------------------#
| IloNumArray | Contains the prices (negative of dual values)                  |
#-------------#----------------------------------------------------------------#
| ms_graph    | The current graph                                              |
#-------------#----------------------------------------------------------------#
| vector<int> | The selected column (a column containing 2 and 4 from a problem|
|             | with 5 nodes looks like that : [0 0 1 0 1])                    |
#-------------#----------------------------------------------------------------#
Once your function has the proper format, call BnP::setGenerator(Generator gen),
with your function pointer as 'gen'.

>> How to use your own branch-selecting methods?
Your methods must have the Generator form: 
(ms_graph, vector<set<int>>, Node)-> pair<int,int> with :
#-------------------#----------------------------------------------------------#
| ms_graph          | The graph in its current state                           |
#-------------------#----------------------------------------------------------#
| vector<set<int>>  | The vector containing every columns as sets              |
#-------------------#----------------------------------------------------------#
| Node              | The _currentNode                                         |                
#-------------------#----------------------------------------------------------#
| pair<int,int>     | Return the selected nodes or (0,0) if nothing has been   |
|                   | found                                                    |
#-------------------#----------------------------------------------------------#

== TUNE THE MS_GRAPH ===========================================================

>> How to choose the node-selecting function for the max set?
You have the choice between 3 functions :
- SN_MODE_MIN_SWAP : Choose the best node to make the intstack update easier.
- SN_MODE_NEAR : Choose the first node with a negative price.
- SN_MODE_SORT : Choose the maximum node by a given node-sorting function that
you can also set.

>> How to use your own node-sorting function?
Your methods must have the NodeSorter form: (ms_Node, ms_Mode)-> bool
It must correspond to the requirements of std::sort().
Then call ms_graph::ms_set_node_sorter() if you have access to the graph
or call BnP::setGraphNodeSorter().

>> How to set the objective value?
Call ms_graph::ms_set_obj().

== TUNE THE CPLEX MASTER SOLVER ================================================

>> How to limit the output of CPLEX on the shell?
Call BnP::setDiscreetMode().
To reverse it call BnP::setNoisyMode().

>> How to limit the CPLEX solver to only 1 thread?
Call BnP::setMonoMode().
To reverse it call BnP::setPolyMode();

== RUN THE SOLVER ==============================================================

>> What are the argument of BnP::run()?
#-----------#------------------------------------------------------------------#
| timelimit | timelimit which trigger a soft timeout (stop at the next         |
|           | checkpoint). It is expressed in seconds. -1 will not apply       |
|           | limit.                                                           |
#-----------#------------------------------------------------------------------#
| log       | Append the final solution to a .csv file.                        | 
#-----------#------------------------------------------------------------------#

>> How to retrieve a solution?
Call BnP::print(OutFormat out)
with out among O_STD, O_SOL, O_DOT.
/!\ O_SOL not implemented yet /!\

== WHAT ELSE IS IN BNP =========================================================

>> Enumerations
#------------------------#-----------------------------------------------------#
| InFormat and OutFormat | Input reading and output writing                    |
#------------------------#-----------------------------------------------------#
| NodeStatus             | Status used by forward                              |            
#------------------------#-----------------------------------------------------#
| Transform              | Transform used to create a node                     |
#------------------------#-----------------------------------------------------#

>> BnP::forward(int u, int v)
Create a child node to *_currentNode through (t,u,v). The new pNode replace 
_currentNode and is added to _nodes. The choice of t in Tranform depends of the 
_currentNode status.
If u and v are not given by the user, they are determined either by taking the
previous values or by calling the selection function. If previous values do 
exist, they are used even if the user does fill values.

NS_CREATED or NS_SOLVED => T_LINK
NS_1C                   => T_MERGE
NS_2C                   => do nothing and put a warning in cout.

>> BnP::backward()
Pop the _currentNode from _nodes and free the allocated memory if it's not the
_incumbent or _rootNode. The _currentNode is the parent node of the previous 
_currentNode.
Do nothing if _currentNode is the root.

== HOW MS_GRAPH FIND SETS ======================================================

For now it creates a searching tree and only evaluate the leaves. Through 
backtrack() and deep_forward(), it searches the first set with a total price 
under the objective.

It was supposed to be stopped a prune branches where the estimated score cannot
is over the objective (Estimation is a lower bound). It has been commented out
as it looks like it lengthen the search more than anything else. It may come 
from an errore in ms_cplt(). 




