-- this file is licenced for use under the GNU General Public Licence version 2

newPackage(
	"Posets",
    	Version => "0.1", 
    	Date => "July 16, 2008",
    	Authors => {
	     {Name => "Sonja Mapes", Email => "mapes@math.columbia.edu", HomePage => "http://www.math.columbia.edu/~mapes/"},
	     {Name => "Gwyn Whieldon", Email => "whieldon@math.cornell.edu", HomePage => "http://www.math.cornell.edu/People/Grads/whieldon.html"},
	     {Name => "Josephine Yu", Email => "jyu@math.mit.edu", HomePage => "http://www-math.mit.edu/~jyu/"}
	     },
    	Headline => "Package for processing posets and order complexes",
    	DebuggingMode => true
    	)
 
export {
     Poset,
     poset,
     DirectedGraph,
     directedGraph,
     adjacencyMatrix,
     allPairsShortestPath,
     transitiveClosure,
--     FullRelationMatrix,
     RelationMatrix,
     compare,
     orderIdeal,
     filter,
     Relations,
     GroundSet,
     Vertices,
     DirectedEdges,
     posetMeet, 
     meetExists, -- do tests for
     posetJoin,
     joinExists,
     isLattice,
     --lcm,
     lcmLattice
     }



Poset = new Type of HashTable

poset = method()
poset(List,List) := (I,C) ->
     new Poset from {
	 symbol GroundSet => I,
	 symbol Relations => C,
     	 symbol RelationMatrix => transitiveClosure(I,C),
	 symbol cache => CacheTable
	 }
    
-- in case you actually have M to begin with    
poset(List,List,Matrix) := (I,C,M) ->
     new Poset from {
	  symbol GroundSet => I,
	  symbol Relations => C,
	  symbol RelationMatrix => M,
	  symbol cache => CacheTable
	  }
     
DirectedGraph = new Type of HashTable
 
directedGraph = method()
directedGraph(List, List) := (V,E) ->
     new DirectedGraph from {
     	  symbol Vertices => V,
     	  symbol DirectedEdges => E,
	  symbol cache => CacheTable
     }

--------------

--inputs: (I,C), I is a List (ground set or vertex set) and C is a List of pairs of elements in I
--     	     OR DirectedGraph OR Poset
--output: a matrix whose rows and columns are indexed by I, 
--     	    where (i,j) entry is infinity (i.e. 1/0.) if (i,j) is not in C
--     	    and 1 otherwise (i.e. tropicalization of the "usual" adjacency matrix)
--caveat: diagonal entries are 0
-- uses:  transitive closure 

adjacencyMatrix = method()
adjacencyMatrix(List,List) := Matrix => (I, C) -> (
     M := mutableMatrix table(#I, #I, (i,j)->1/0.);
     ind := hashTable( apply(I, i-> i=> position(I,j-> j==i) )  ); --indices 
     scan(C, e -> M_(ind#(e#0), ind#(e#1))= 1);
     scan(numrows M, i-> M_(i,i) = 0);
     matrix M      
     )
adjacencyMatrix(DirectedGraph) := Matrix => (G) -> adjacencyMatrix(G.Vertices,G.DirectedEdges)
adjacencyMatrix(Poset) := Matrix => (P) -> adjacencyMatrix(P.GroundSet,P.Relations)

--input: adjacency matrix of a directed graph
--output: a matrix whose (i,j) entry is the length of the shortest path from i to j
--algorithm: Floyd–Warshall algorithm for all pairs shortest path
allPairsShortestPath = method()
allPairsShortestPath(Matrix) := Matrix => (A) -> (
     D := mutableMatrix(A);
     n := numrows D;
     scan(n, k->
	       table(n,n,(i,j)-> D_(i,j) = min(D_(i,j), D_(i,k)+D_(k,j)))
	  );
     matrix D
     )
allPairsShortestPath(DirectedGraph) := Matrix => (G)-> allPairsShortestPath(adjacencyMatrix(G))



-- input: a poset, and an element A from I
-- output:  the index of A in the ground set of P
-- usage: compare, orderIdeal 
indexElement := (P,A) -> (
      sum apply(#P.GroundSet, i-> if P.GroundSet#i == A then i else 0))

-- input:  a list, potentially with nulls
-- output:  a list w/out nulls
-- usage:  orderIdeal, filter
nonnull :=(L) -> (
     select(L, i-> i =!= null))


--------------------------------------------------
--Transitive Closure and Element Inclusion
--------------------------------------------------

--input: (I,C).  I=List of vertices.  C=List of pairs (edges)
--output: matrix where 1 in (i,j) position where i <= j, 0 otherwise
--uses: poset

transitiveClosure = method()
transitiveClosure(List,List) := List => (I,C)-> (
     A := adjacencyMatrix(I,C);
     D := mutableMatrix allPairsShortestPath(A);
     scan(numrows D, i-> D_(i,i) = 0);
     table(numrows D, numrows D, (i,j)->(
	  if D_(i,j) ==1/0. then D_(i,j) = 0 else D_(i,j) = 1;	  
	  )
	  );
     matrix D
     )



--input: A poset with any type of relation C (not necessarily minimal, maximal, etc.)
--output:  The transitive closure of relations in C in our poset

--fullPosetRelation:= (P) -> (
--     M:=P.RelationMatrix;
--     L = toList sum apply(numrows(M), i-> set(nonnull(apply(numrows(M), 
--	       j-> if (M_j)_i=!=0 and i=!=j then (P.GroundSet#i,I#j)))))
--     ) -- outputs the wrong data

--input: A poset P with any type of relation C (minimal, maximal, etc.)
--output:  The poset P' on the same ground set with the transitive closure of C

--fullPoset:= (P) -> (
--     L = poset(P.GroundSet,fullPosetRelation(P)) 
--     )

-- input:  A poset, and two elements A and B from I
-- output: true if A<= B, false else
compare = method()
compare(Poset, Thing, Thing) := Boolean => (P,A,B) -> (
     Aindex:=indexElement(P,A);
     Bindex:=indexElement(P,B);
          if P.RelationMatrix_Bindex_Aindex==0 then false else true
     )



--------------------------------------------------
--Covering Relations
--------------------------------------------------

testcover=(P,A,B) -> (
     L:=poset(P.GroundSet,fullPosetRelation(P));
     k:=#L.GroundSet-2; 
         
     if sum(nonnull(apply(k, i-> if compare(L,A,(toList(set(L.GroundSet)-{A,B}))_i)==true and
	       compare(L,(toList(set(L.GroundSet)-{A,B}))_i,B)==true
	       then 1)))=!=0 then C=C+set{(A,B)};
     C
)		

--input: A poset with any type of relation C (minimal, maximal, etc.)
--output: The minimal relations defining our poset

coveringRelations:=(P) -> (
     C=set{};
     apply(#P.CRelations,i->testcover(P,P.CRelations#i#0,P.CRelations#i#1));
     toList(set(P.CRelations)-C)
     )

--input: A poset with any type of relation C (minimal, maximal, etc.)
--output:  A new poset P with the minimal relations

coveringRelationsPoset:=(P) -> (
     L=poset(P.GroundSet,coveringRelations(P))
     )

--------------------------------------------------
--Minimal Element Construction
--------------------------------------------------

minimalElementIndex:=(P)-> (
     M:=P.RelationMatrix;
     nonnull(apply(numcols(M), k-> if (apply(numcols(M), j-> (sum((apply(numrows(M),i-> (transpose(M))_i))))_j))#k==1 then k))
     )

minimalElements:=(P) -> (
     L:=minimalElementIndex(P);
     apply(#L,i-> P.GroundSet#(L#i))
     )

PosetMinusMins:=(P)-> (
     L:=minimalElements(P);
     K:=fullPoset(P);
     N:=set{};
     S:=apply(#L, j-> apply(#K.CRelations,i->(K.CRelations#i)#0===L#j));
     E:=sum set nonnull(apply(#K.CRelations,l->if member(true,set apply(#L,k->S#k#l)) then N=N+set{K.CRelations#l}));
     C:=toList (set(K.CRelations)-N);
     I:=toList (set(K.GroundSet)-set(L));
     poset(I,C)
     )



--------------------------------------------------
--Order and Filter Ideals
--------------------------------------------------

-- input: a poset, and an element from I
-- output: the order ideal of a, i.e. all elements in the poset that are >= a
-- usage:

orderIdeal= method()
orderIdeal(Poset, Thing) := (P, a) -> (
     M:=P.RelationMatrix;
     aindex := indexElement (P,a);
     GreaterThana:= entries((transpose(M))_aindex);
     nonnull(apply(#GreaterThana, i-> if GreaterThana_i == 1 then P.GroundSet#i)) 
     )	
	   


-- input: a poset, and an element from I
-- output:  the filter of a, i.e. all elements in the poset that are <= a
-- usage:
filter = method()
filter(Poset, Thing) := (P,a) -> (
     M:=P.RelationMatrix;
     aindex := indexElement (P,a);
     LessThana:= entries M_aindex;
     nonnull(apply(#LessThana, i-> if LessThana_i == 1 then P.GroundSet#i))
     )



----------------------------------------------------
--Joins, Meets, Lattices and Atoms
----------------------------------------------------
-- inputs: P, poset, and two elements of P.GroundSet
-- outputs:  the element of P.GroundSet that is the join of these, or "not comparable" or "not unique" if those situations occur
-- usage:  joinExists used in isLattice

posetJoin = method()     
posetJoin(Poset,Thing,Thing) := (P,a,b)  -> (
     OIa := orderIdeal(P,a);     
     OIb := orderIdeal(P,b);
     upperBounds := toList (set(OIa)*set(OIb));
     if upperBounds == {} then (error "your elements do not share any upper bounds") else (M := P.RelationMatrix;
     	  heightUpperBounds := flatten apply(upperBounds, element-> sum entries M_{indexElement(P,element)});
     	  if #(select(heightUpperBounds, i-> i== min heightUpperBounds)) > 1 then error "join does not exist, least upper bound not unique" 
	  else(upperBounds_{position (heightUpperBounds, l -> l == min heightUpperBounds)})
     ))


joinExists = method()
joinExists(Poset,Thing,Thing) := (P,a,b) -> (
     OIa := orderIdeal(P,a);     
     OIb := orderIdeal(P,b);
     upperBounds := toList (set(OIa)*set(OIb));
     if upperBounds == {} then false else (M := P.RelationMatrix;
     	  heightUpperBounds := flatten apply(upperBounds, element-> sum entries M_{indexElement(P,element)});
     	  if #(select(heightUpperBounds, i-> i== min heightUpperBounds)) > 1 then false
	  else true
     ))


--inputs:  P a poset, and 2 elements of P.GroundSet
--outputs:  the element in P.GroundSet that is the meet of these, or "not comparable" or "not unique"
-- usage:  MeetExits used in isLattice
posetMeet = method()
posetMeet(Poset,Thing,Thing) := (P,a,b) ->(
     Fa:= filter(P,a);
     Fb:= filter(P,b);
     lowerBounds:= toList (set(Fa)*set(Fb));
     if lowerBounds == {} then (error "your elements do not share any lower bounds") else (
	  M := P.RelationMatrix;
     	  heightLowerBounds := flatten apply(lowerBounds, element-> sum entries M_{indexElement(P,element)});
     	  if #(select(heightLowerBounds, i-> i== max heightLowerBounds)) > 1 then error "meet does not exist, greatest lower bound not unique" 
	  else(lowerBounds_{position (heightLowerBounds, l -> l == max heightLowerBounds)})) 
     )

meetExists = method()
meetExists(Poset, Thing, Thing) := (P,a,b) -> (
     Fa:= filter(P,a);
     Fb:= filter(P,b);
     lowerBounds:= toList (set(Fa)*set(Fb));
     if lowerBounds == {} then false else (
	  M := P.RelationMatrix;
     	  heightLowerBounds := flatten apply(lowerBounds, element-> sum entries M_{indexElement(P,element)});
     	  if #(select(heightLowerBounds, i-> i== max heightLowerBounds)) > 1 then false
	  else true ) 
     )


--inputs: a poset P
--output:  boolean value for whether or not it is a lattice
--usage:
isLattice = method()
isLattice(Poset) := (P) -> (
    checkJoins := unique flatten flatten apply(P.GroundSet, elt -> 
	                apply (P.GroundSet, elt2-> joinExists(P,elt, elt2)));
    checkMeets :=  unique flatten flatten apply(P.GroundSet, elt -> 
	                apply (P.GroundSet, elt2-> meetExists(P,elt, elt2) ));
    if member(false, set (flatten{checkJoins,checkMeets})) === true then false else true 
     )


-----------------------------------------------
-- LCM lattices
-----------------------------------------------
--input: a set of monomials
-- output: the lcm of those monomials
--lcm = (L) -> (
--    flatten entries gens intersect apply(L, i-> ideal (i))
--    )

-- input:  generators of a monomial ideal
-- output: lcm lattice of that monomial ideal, without the minimal element
-- potential problem:  subsets dies when a set is too big (> 18)

lcmLattice = method()     
lcmLattice(Ideal) := Poset => (I) -> (
	   L := flatten entries gens I;
	   subsetsL := flatten apply(#L, i-> subsets (L,i+1));
	   Ground := unique flatten apply (subsetsL, r-> lcm(r));
	   Rels := nonnull unique flatten apply (Ground, r-> apply(Ground, s-> if s%r == 0 then (r,s)));
	   RelsMatrix :=  matrix apply (Ground, r-> apply(Ground, s-> if s%r == 0 then 1 else 0));
	   P = poset (Ground, Rels, RelsMatrix);
	   P)

lcmLattice(MonomialIdeal) := Poset => (M) -> (
     	   L := flatten entries gens M;
	   subsetsL := flatten apply(#L, i-> subsets (L,i+1));
	   Ground := unique flatten apply (subsetsL, r-> lcm(r));
	   Rels := nonnull unique flatten apply (Ground, r-> apply(Ground, s-> if s%r == 0 then (r,s)));
	   RelsMatrix :=  matrix apply (Ground, r-> apply(Ground, s-> if s%r == 0 then 1 else 0));
	   P = poset (Ground, Rels, RelsMatrix);
	   P)
----------------------------------
-- DOCUMENTATION
----------------------------------

beginDocumentation()


---------
-- front page
---------
doc ///
     Key
     	  Posets
     Headline
          A package for working with posets. 
     Description
          Text
	       {\em Posets} package defines Poset as a new data type and provides 
	       routines which use or produce posets.   A poset or a partially ordered set is a set together with a binary relation satisfying reflexivity, antisymmetry, and transitivity.
///
	 
----------
-- types
----------	       
doc ///
     Key
     	  Poset
     Headline
     	  a class for partially ordered sets (posets)
     Description
     	  Text
     	       This class is a type of HashTable which represents finite posets.  It consists of a ground 
	       set, a set of sequences (a,b) representing relationships where a is less than or 
	       equal to b, a matrix encoding these relations. 
	  Example
	       G = {a,b,c,d,e}; -- the ground set
	       R = {(a,b),(b,c),(a,c),(a,d),(d,e)}; --a set of sequences representing relations that "generate" all other relations
	       P = poset (G,R) -- the matrix encoding these relations is computed by calling this function
	       
     SeeAlso
     	  poset
	  GroundSet
	  Relations
	  RelationMatrix	   
///     


doc ///
     Key
     	  DirectedGraph
     Headline 
     	  a class for directed graphs
     Description
     	  Text
	       This class is a type of HashTable which represents directed graphs.  It consists of a vertex 
	       set and a set of sequences (a,b) representing a directed edge from a to b. 
	  Example
	       V = {a,b,c,d,e}; -- the vertex set
	       E = {(a,b),(b,c),(a,c),(a,d),(d,e)}; -- directed edges; (a,b) means a directed edge from a to b.
	       G = directedGraph(V,E)
	       allPairsShortestPath(G)
     SeeAlso
     	  directedGraph
	  Vertices
	  DirectedEdges
	  adjacencyMatrix
	  allPairsShortestPath
/// 


-------------
-- constructors
-------------

doc ///
     Key
     	  poset
	  (poset, List, List)
	  (poset, List, List, Matrix)
     Headline
     	  creating a poset
     Usage
     	   P = poset (G,R)
	   P = poset (G,R,M)
     Inputs
     	  G : List 
	       elements in the ground set of P
	  R : List 
	       sequences (a,b) which indicate that a is less than or equal to b
	  M : Matrix
	       with entries ij equal to one when the jth element in G is less than or equal to the ith element in G
     Outputs
     	  P : Poset
	       a poset consisting of the elements in G, with the order relations determined by R and or M
     Description
      	   Text
	   	This function allows one to create a poset by defining the set and giving the order relations between the elements in the set.
	   Example
	   	G = {a,b,c,d};
		R = {(a,b), (a,c), (c,d)};
		P = poset (G,R)
	   Text
	   	Sometimes in finding "enough" relations one inadverdently finds all of the relations in the poset.  In this case, the user
		can bypass having the function poset calculate the transitive closure of the relations by entering the matrix encoding the
		relations in when calling the poset function.
	   Example
	   	S = QQ[x,y,z];
		G = {x^2, x*y, z^2, x^2*y*z, x*y*z^3, x^2*y^2*z^3};
	        R = select((flatten apply (G, g-> apply (G, h-> if h % g == 0 then (g,h)))), i -> i =!= null) -- finds all pairs where g divides h
		M = matrix apply (G, g-> apply (G, h-> if h % g == 0 then 1 else 0)) 
		P = poset(G,R,M) 
		
     SeeAlso
     	  Poset
	  GroundSet
	  Relations
	  RelationMatrix
///
 
doc ///
     Key
     	  directedGraph
	  (directedGraph, List, List)
     Headline
     	  creating a directed graph
     Usage
     	  G = directedGraph(V,E)
     Inputs
     	  V : List 
	       elements in the vertices of G
	  E : List 
	       sequences (a,b) which represents a directed edge from a to b
     Outputs
     	  G : DirectedGraph 
	       a directed graph on vertex set V with directed edges E
     Description
      	   Text
	   	This function allows one to create a directed graph by defining the vertices and directed edges between them.
	   Example
	   	V = {a,b,c,d};
		E = {(a,b), (a,c), (c,d)};
		G = directedGraph (V,E)	
     SeeAlso
     	  DirectedGraph
	  Vertices
	  DirectedEdges
     	  adjacencyMatrix     
/// 
 
-----------
-- finding relations
-----------

doc ///
     Key 
     	  adjacencyMatrix
	  (adjacencyMatrix, DirectedGraph)
	  (adjacencyMatrix, Poset)
	  (adjacencyMatrix, List, List)
     Headline
	  returns adjacency matrix of a directed graph
     Usage
	  M = adjacencyMatrix(G)
	  M = adjacencyMatrix(P)
	  M = adjacencyMatrix(I,C)
     Inputs
	  G : DirectedGraph
	  P : Poset
       	  I : List
	  C : List
     Outputs
      	  M : Matrix 
	       whose rows and columns are indexed by G.Vertices 
	       or P.GroundSet or I. The (i,j) entry of M  is 1 if (i,j) is 
	       in G.DirectedEdges or P.Relations or C, 
	       and infinity otherwise.  Diagonal entries are 0. 
     Description
	  Example
     	       I = {a,b,c,d,e};
	       C = {(a,b),(b,c),(a,c),(a,d),(d,e)};
	       G = directedGraph(I,C);
	       P = poset(I,C);
	       adjacencyMatrix(G)
	       adjacencyMatrix(P)
	       adjacencyMatrix(I,C)	       	       
     Caveat
     	  Diagonal entries are 0.  Output matrix is over RR.
     SeeAlso
     	  DirectedGraph
     	  allPairsShortestPath    
/// 
 
doc ///
     Key 
     	  allPairsShortestPath
	  (allPairsShortestPath, DirectedGraph)
	  (allPairsShortestPath, Matrix)
     Headline
     	  computes lengths of shortest paths between all pairs in a directed graph
     Usage
     	  D = allPairsShortestPath(G)
	  D = allPairsShortestPath(A)
     Inputs
     	  G : DirectedGraph 
	  A : Matrix
	      adjacency matrix of a directed graph, whose (i,j) entry is the length of the directed edge from i to j and is infinity (1/0.) if there is no such directed edge.
     Outputs
     	  D : Matrix 
	       whose (i,j) entry is the length of the shortest path from i to j and is infinity if there is no directed path from i to j.
     Description
      	   Text
	   	This function uses the Floyd-Warshall algorithm to compute the lengths of shortest paths between all pairs of vertices in a directed graph.
	   Example
	       G = directedGraph( {a,b,c,d,e}, {(a,b),(b,c),(a,c),(a,d),(d,e)});
	       allPairsShortestPath(G)
	       allPairsShortestPath(adjacencyMatrix(G))
	       A = matrix({{0,1,3,5},{1/0.,0,1,3},{1/0.,1/0.,0,1},{2,1/0.,1/0.,0}})
	       allPairsShortestPath(A)
     Caveat
     	  Assume there is no negative cycles.  Output matrix is over RR.
     SeeAlso
     	  DirectedGraph
     	  adjacencyMatrix   
/// 
  

doc ///
     Key
     	  transitiveClosure
	  (transitiveClosure, List, List)
     Headline
     	  computes the transitive closure of a given set of relations.  
     Usage
     	  M = transitiveClosure(I,R)
     Inputs
     	  I : List
	       ground set 
	  R : List
	       of pairs of elements in the ground set I.
     Outputs
     	  M : Matrix 
	       whose (i,j) entry is 1 if (i,j) is a relation and 0 otherwise.
     Description
      	   Text
	   	This function uses allPairsShortestPath method is used by the poset constructor to compute RelationMatrix from Relations in a Poset.
	   Example
	       I = {a,b,c,d,e}; -- the ground set
	       R = {(a,b),(b,c),(a,c),(a,d),(d,e)}; -- relations
     	       transitiveClosure(I,R)
     Caveat
     	  Output matrix is over RR.
     SeeAlso
     	  Poset
	  poset
	  Relations
	  RelationMatrix
     	  allPairsShortestPath
/// 
 
 
---------
-- poset routines
---------
 
doc ///
     Key
     	  compare
	  (compare, Poset,Thing,Thing)
     Headline
	  returns boolean value for whether an element is less than another 
     Usage
	  compare(P,a,b)
     Inputs
     	  P : Poset
	  a : Thing
	  b : Thing
	       a and b are elements of P
     Outputs
      	   true : Boolean
	   	if a is less than or equal to b
	   false : Boolean
	   	otherwise
     Description
     	  Text
	       This function simply looks at two elements in P, and determines what if any relation exists between them.
	  Example
	       P = poset ({a,b,c}, {(a,b), (a,c)});
	       compare (P,a,c)
	       compare (P,c,a)
	       compare (P,b,c)
	       
	       
     Caveat
     	  Note that if two elements are not comparable, compare still returns a value of false. 
	  
/// 
  
 
doc ///
     Key
     	  filter
	  (filter,Poset,Thing)
     Headline
     	  returns a principal filter generated by the given element
     Usage 
     	  F = filter (P,a)
     Inputs
     	  P : Poset
	  a : Thing
	       a is an element of P
     Outputs
     	  F : List
	       a list of all elements in P that are less than or equal to a
     Description
     	  Example
	       G = {a,b,c,d};
	       R = {(a,b), (a,c), (c,d)};
	       P = poset (G,R);
	       F = filter (P,d)
     SeeAlso
     	  orderIdeal	       
/// 

doc ///
     Key
     	  orderIdeal
	  (orderIdeal, Poset, Thing)
     Headline
     	  returns a principal order ideal generated by the given element
     Usage 
     	  O = orderIdeal (P,a)
     Inputs
     	  P : Poset
	  a : Thing
	       a is an element of P
     Outputs
     	  O : List
	       a list of all elements in P that are greater than or equal to a
     Description
     	  Example
	       G = {a,b,c,d};
	       R = {(a,b), (a,c), (c,d)};
	       P = poset (G,R);
	       O = filter (P,c)
     SeeAlso
     	  filter    
/// 
 
 
doc ///
     Key
     	  posetMeet
	  (posetMeet, Poset, Thing, Thing)
     Headline
     	  returns the meet of two elements
     Usage
          m = posetMeet (P,a,b)
     Inputs
     	  P : Poset
	  a : Thing
	  b : Thing
	       a and b are in P
     Outputs
     	  m : Thing
	       m is the meet of a and b in P
     Description 
     	  Text
	       This function returns the greatest element that is less than both a and b in P.
	  Example
	       P = poset ({a,b,c,d,e,f}, {(a,d),(d,f),(b,d),(b,e),(c,e),(e,f)});
	       posetMeet (P, d,e)
     Caveat
      	   This function returns an error if the meet does not exist.  To check if the meet exists use meetExists.
     SeeAlso
      	   meetExists
	   posetJoin
	   joinExists
/// 



doc ///
     Key
     	  posetJoin
	  (posetJoin, Poset, Thing, Thing)
     Headline
     	  returns the join of two elements
     Usage
          j = posetJoin (P,a,b)
     Inputs
     	  P : Poset
	  a : Thing
	  b : Thing
	       a and b are in P
     Outputs
     	  j : Thing
	       j is the join of a and b in P
     Description 
     	  Text
	       This function returns the least element that is greater than both a and b in P.
	  Example
	       P = poset ({a,b,c,d,e,f}, {(a,d),(d,f),(b,d),(b,e),(c,e),(e,f)});
	       posetJoin (P, a,b)
     Caveat
      	   This function returns an error if the join does not exist.  To check if the join exists use joinExists.
     SeeAlso
      	   joinExists
	   posetMeet
	   meetExists
/// 
   
   
doc ///
     Key
     	  meetExists
	  (meetExists, Poset, Thing, Thing)
     Headline
     	  determines if the meet exists
     Usage
          meetExists (P,a,b)
     Inputs
     	  P : Poset
	  a : Thing
	  b : Thing
	       a and b are in P
     Outputs
     	  true : Boolean
	       the meet of a and b in P exists
	  false : Boolean
	       otherwise
     Description 
     	  Text
	       This function determines if greatest element that is less than both a and b in P exists.
	  Example
	       P = poset ({a,b,c,d,e,f}, {(a,d),(d,f),(b,d),(b,e),(c,e),(e,f)});
	       meetExists (P,d,e)
     	       meetExists(P,a,b)
     Caveat
      	  If the meet exists, to find it use posetMeet.
     SeeAlso
      	   posetMeet
	   posetJoin
	   joinExists    
/// 

doc ///
     Key
     	  joinExists
	  (joinExists, Poset, Thing,Thing)
     Headline
     	  determines if the join exists
     Usage
          joinExists (P,a,b)
     Inputs
     	  P : Poset
	  a : Thing
	  b : Thing
	       a and b are in P
     Outputs
     	  true : Boolean
	       the join of a and b in P exists
	  false : Boolean
	       otherwise
     Description 
     	  Text
	       This function determines if least element that is greater than both a and b in P exists.
	  Example
	       P = poset ({a,b,c,d,e}, {(a,d),(b,d),(b,e),(c,e)});
	       joinExists (P,d,e)
     	       joinExists(P,a,b)
     Caveat
      	  If the join exists, to find it use posetJoin.
     SeeAlso
      	   posetJoin
	   posetMeet
	   meetExists  
/// 
 
doc ///
     Key
     	  isLattice
	  (isLattice,Poset)
     Headline
     	  determines if a poset is a lattice
     Usage
     	  isLattice (P)
     Inputs
     	  P: Poset
     Outputs
     	  true : Boolean
	       if P is a lattice
	  false : Boolean
	       otherwise
     Description 
     	  Text
	       This function examines the entire poset to determine whether or not every pair of elements has both a meet and a join.
	  Example
	       P = poset ({a,b,c,d,e,f}, {(a,d),(d,f),(b,d),(b,e),(c,e),(e,f)});
	       isLattice (P)
	  Text
	       And by adding an element to the above example, we can create a poset which is a lattice.     
	  Example
	       P = poset ({a,b,c,d,e,f,x}, {(a,d),(d,f),(b,d),(b,e),(c,e),(e,f), (x,a), (x,b), (x,c)});   
    	       isLattice (P)	   
    
/// 

     
---------
-- LCM lattices
---------

doc ///
     Key
     	  lcmLattice
	  (lcmLattice,Ideal)
	  (lcmLattice, MonomialIdeal)
     Headline
     	  returns the LCM lattice of an ideal
     Usage
     	  L = lcmLattice (I)
	  L = lcmLattice (M)
     Inputs 
     	  I : Ideal
	  M : MonomialIdeal
     Outputs
     	  L : Poset
     Description
     	  Text
	       This command allows for fast computation of LCM lattices, which are particularly useful in the study of resolutions of monomial ideals.
	       Specifically the LCM lattice is the set of all lcms of subsets of the generators of the ideal, partially ordered by divisability.  
	  Example
	       S = QQ[a,b,c,d];
	       I = ideal (b^2-a*d, a*d-b*c, c^2-b*d);
	       M = monomialIdeal (b^2, b*c, c^2);
	       L = lcmLattice (I);
	       L.GroundSet
	       L.RelationMatrix
	       LM = lcmLattice (M);
	       LM.GroundSet
	       LM.RelationMatrix
     Caveat           
     	  Note, that at present this command does not efficiently handle ideals with large numbers of generators.  This is a problem that should be
	  fixed by the next release of this package.	 
/// 

doc ///
     Key 
     	  GroundSet
     Headline
     	  underlying set of a poset
     Usage
     	  G = P.GroundSet
     Inputs
     	  P : Poset
     Outputs
     	  G : List
	       list of elements in P without the data of how they are partially ordered
     Description
     	  Text
     	       Since any poset is in fact a HashTable this symbol denotes the data in the HashTable consisting of the elements of the set.
	  Example
	       S = QQ[a,b,c,d];
	       M = monomialIdeal (b^2, b*c, c^2);
	       L = lcmLattice (M);
	       L.GroundSet
     SeeAlso
     	  RelationMatrix
	  Relations
	  poset
	  Poset
///
 
 
doc ///
     Key 
     	  RelationMatrix
     Headline
     	  the matrix expressing all of the relations between elements in a Poset
     Usage
     	  M = P.RelationMatrix
     Inputs
     	  P : Poset
     Outputs
     	  M : Matrix
	       where the ijth entry indicates whether or not the jth element in the set is less than or equal to the ith element
     Description
     	  Text
     	       Since any poset is in fact a HashTable this symbol denotes the data in the HashTable containing all possible relations between elements.
	  Example
	       S = QQ[a,b,c,d];
	       M = monomialIdeal (b^2, b*c, c^2);
	       L = lcmLattice (M);
	       L.GroundSet
	       L.RelationMatrix
     SeeAlso
     	  GroundSet
	  Relations
	  poset
	  Poset
/// 
 
doc ///
     Key 
     	  Relations
     Headline
     	  a set of relations in the poset that generates all other relations
     Usage
     	  R = P.Relations
     Inputs
     	  P : Poset
     Outputs
     	  R : List
	       list of pairs of elements in P where (a,b) means a is less than or equal to b
     Description
     	  Text
     	       Since any poset is in fact a HashTable this symbol denotes the data in the HashTable which will describe all of the relations
	       between elements in P.
	  Example
	       P = poset ({a,b,c,d,e,f},{(a,d),(d,f),(b,d),(b,e),(c,e),(e,f)});
	       P.Relations --note there is not a one to one correspondence between these and entries in the RelationMatrix
	       P.RelationMatrix
     Caveat
     	  Typically, the relations are the data entered in by the user and by no means account for every relation.  The RelationMatrix is usually computed
	  from this set and thus is a better descriptor of the relations in P.
     SeeAlso
     	  RelationMatrix
	  GroundSet
	  poset
	  Poset
///

doc ///
     Key 
     	  Vertices
     Headline
     	  the set of vertices a directed graph
     Usage
     	  V = G.Vertices
     Inputs
     	  G : DirectedGraph
     Outputs
     	  V : List
	       list of vertices of the directed graph G
     Description
     	  Text
     	       Since any DirectedGraph is in fact a HashTable this symbol denotes the data in the HashTable consisting of the vertices of the directed graph.
	  Example
	       G = directedGraph({a,b,c,d,e},  {(a,b),(b,c),(a,c),(a,d),(d,e)});
	       V = G.Vertices -- the vertex set
     SeeAlso
     	  DirectedGraph
	  directedGraph
	  DirectedEdges
///

doc ///
     Key 
     	  DirectedEdges
     Headline
     	  the set of directed edges of a directed graph
     Usage
     	  E = G.DirectedEdges
     Inputs
     	  G : DirectedGraph
     Outputs
     	  E : List
	       list of directed edges of the directed graph G
     Description
     	  Text
     	       Since any DirectedGraph is in fact a HashTable this symbol denotes the data in the HashTable consisting of the directed edges of the directed graph.
	  Example
	       G = directedGraph({a,b,c,d,e},  {(a,b),(b,c),(a,c),(a,d),(d,e)});
	       E = G.DirectedEdges -- the vertex set
     SeeAlso
     	  DirectedGraph
	  directedGraph
	  Vertices
///
 
 
 
---------------------------------
--Tests
---------------------------------

-- a lattice, B_3
TEST ///
I ={a,b,c,d,e,f,g,h};
C ={(a,b),(a,c),(a,d),(b,e),(b,f),(c,e),(c,g),(d,f),(d,g),(e,h),(f,h),(g,h)};
P=poset(I,C);
M = matrix {{1,1,1,1,1,1,1,1},
	    {0,1,0,0,1,1,0,1},
	    {0,0,1,0,1,0,1,1},
	    {0,0,0,1,0,1,1,1},
	    {0,0,0,0,1,0,0,1},
	    {0,0,0,0,0,1,0,1},
	    {0,0,0,0,0,0,1,1},
	    {0,0,0,0,0,0,0,1}};
assert (entries P.RelationMatrix == entries M)
--G=directedGraph(I,C)
--A=adjacencyMatrix(I,C) -- not exported
--allPairsShortestPath(A) -- not exported
--adjacencyMatrix(G) -- not exported
--adjacencyMatrix(P) -- not exported
--transitiveClosure(I,C) 
assert (posetJoin(P,a,b) == {b})
assert (posetJoin(P,b,d) == {f})
assert (posetMeet(P,a,b) == {a})
assert (posetMeet(P,f,g) == {d})
assert (orderIdeal(P,a) == {a,b,c,d,e,f,g,h})
assert (orderIdeal(P,b) == {b,e,f,h})
assert (filter(P,a) == {a})
assert (filter(P,g) == {a,c,d,g})
assert (isLattice(P))
///


-- two equivllaent non lattices with different initial data
TEST ///
I1={a,b,c,d,e,f};
C1={(a,c),(a,d),(b,c),(b,d),(c,e),(d,e),(e,f)};
P1=poset(I1,C1);
--G1 = directedGraph(I1,C1)

--Poset P1 with additional relations (a,e) and (a,f) added
I2={a,b,c,d,e,f};
C2={(a,c),(a,d),(b,c),(b,d),(c,e),(d,e),(a,e),(a,f),(e,f)};
P2=poset(I2,C2);
--G2=directedGraph(I2,C2) -- adds the non minimal edges

assert (P1.RelationMatrix == P2.RelationMatrix)    
assert (filter(P1,b) == {b})
assert (filter(P1,c) == {a,b,c})
assert (orderIdeal (P1,b) == {b,c,d,e,f})
assert (isLattice (P1) == false)
-- do joins and meets
///

-- do an LCM lattice
-- do order ideal

TEST ///
V = {a,b,c,d,e};
E = {(a,b),(a,d),(b,c),(b,e),(c,e),(d,e)};
G = directedGraph(V,E);
A = adjacencyMatrix(G);
D = allPairsShortestPath(G);
T = transitiveClosure(V,E);

assert(A_(0,1)==1 and A_(0,3)==1 and A_(1,2)==1 and A_(1,4)==1 and A_(2,4)==1 and A_(3,4)==1)
assert(D_(0,4)==2 and D_(4,0)==1/0. and D_(3,3)==0)
assert(T== promote(matrix {{1, 1, 1, 1, 1}, {0, 1, 1, 0, 1}, {0, 0, 1, 0, 1}, {0, 0, 0, 1, 1}, {0, 0, 0, 0, 1}}, RR))

D1 =allPairsShortestPath(matrix({{0,1,1/0.,4},{1/0.,0,1/0.,2},{1,1/0.,0,1/0.},{1/0.,1/0.,1,0}})); -- digraph with a cycle
assert(D1 ==  promote(matrix {{0, 1, 4, 3}, {4, 0, 3, 2}, {1, 2, 0, 4}, {2, 3, 1, 0}}, RR))

///
