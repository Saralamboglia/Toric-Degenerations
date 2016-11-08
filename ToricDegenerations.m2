newPackage(
    	"ToricDegenerations",
	Version => "0.1",
	Date => "October 2016",
	Authors => {
                   {Name => "Sara Lamboglia", Email => "S.Lamboglia@warwick.ac.uk", HomePage=>""}
		},
	Headline => "Package to compute toric degenerations",
	Configuration => {
		"path" => "",
		"fig2devpath" => "",
		"keepfiles" => false,
		"cachePolyhedralOutput" => true,
		"minConvention" => true
	},
        PackageExports => {"FourierMotzkin","gfanInterface","Binomials","Tropical"},
	DebuggingMode => true
)
export {
"computeWeightVectors",
"getNonPrimeIdeals",
"getNonPrimeVectors",
"getPrimeVectors",
"getNonBinomialIdeals", 
"getNonPrimeCones", 
"getPrimeCones",
"findToricDegenerations",
"findInitialIdeals",
"minConvention",
"mapMaximalCones",
"findMissingBinomials",
"isBinomialIdeal"
}

-------------------------------------------------------------------------------------------------------------
--CODE
-------------------------------------------------------------------------------------------------------------





--Computes one vector in the relative interior of a cone of the form r_1+....+r_d(ausiliary function)
getVector = (mycone,rays) -> (
    sum(apply(mycone, idx->rays_idx))
);
--input: one maximal cone (the list of its rays) and the list of rays
--output: the weight vector for that cone


--Computes  one weight vector for each cone, available MIN or MAX convention
computeWeightVectors=method(Options => {
	minConvention=>true});
computeWeightVectors (List,List,List):=opt->(ConesT,RaysT,LinealitySpaceT)->(
    if opt.minConvention==false then(
     apply(ConesT, mycone->(getVector(mycone,RaysT)-(min(getVector(mycone,RaysT))-1)*sum(LinealitySpaceT))))
     else ( apply(ConesT, mycone->(getVector(mycone,RaysT)-(max(getVector(mycone,RaysT))+1)*sum(LinealitySpaceT))))
     );
--input: maximal cones,rays and linealityspace of a tropical variety , these have to be lists,use the option minConvetion=>false if you are
--using max convention


computeWeightVectors (TropicalCycle):=opt->(T)->(
    RaysT:=getRays T;
    LinealitySpaceT:=getLinealitySpace T;
    ConesT:=getMaximalCones T;
    if opt.minConvention==false then(
     apply(ConesT, mycone->(getVector(mycone,RaysT)-(min(getVector(mycone,RaysT))-1)*sum(LinealitySpaceT))))
     else ( apply(ConesT, mycone->(getVector(mycone,RaysT)-(max(getVector(mycone,RaysT))+1)*sum(LinealitySpaceT))))
     );


--output: a list of weight vectors,one for each maximal cone;in case it's min they are all negative otherwise all positive;

--checks of an ideal is binomial
isBinomialIdeal = ideall -> (
    genss := flatten entries mingens ideall;
    all(genss, poly->(#terms poly) == 2)
);






--Computes if the variety has good cones or not, bad vectors, bad ideals, good vectors, good ideals,badIdeals,badCones and goodCones
findToricDegenerations=method(Options => {
	minConvention=>true});

findToricDegenerations(List,Ideal):=opt->(L,I)->(
tableWithOutput:=new  MutableHashTable;
NonPrimeIdeals:={};
ToricIdeals:={};
PrimeVectors:={};
NonPrimeVectors:={};
NonBinomialIdeals:={};
PrimeCones:={};
NonPrimeCones:={};
initialId:=ideal();
i:=0;
    while (i<#L)  do (
    if opt.minConvention==true then(
    initialId = initialIdeal(-L_i, I))
    else (initialId = initialIdeal(L_i, I));
    --print L_i;
    if  isBinomialIdeal(initialId)==true 
    then (
	if binomialIsPrime(initialId)==true
	then( 
    	    PrimeVectors=append(PrimeVectors,L_i);
	    PrimeCones=append(PrimeCones,i);
	    ToricIdeals=append(ToricIdeals,initialId);
	    print "Found One Good";
	    )
     	else ( 
	    NonPrimeIdeals=append(NonPrimeIdeals,initialId);
            NonPrimeVectors=append(NonPrimeVectors,L_i);
            NonPrimeCones=append(NonPrimeCones,i)
	    );
    ) else (
       NonBinomialIdeals=append(NonBinomialIdeals,initialId);
       print "Attention there is a non binomial initial ideal"   
           )
   ;i=i+1;
);
--(badVectors,badIdeals,goodVectors,goodIdeals,veryBadIdeals)
tableWithOutput#"ToricIdeals" = ToricIdeals;
tableWithOutput#"NonPrimeIdeals" = NonPrimeIdeals;
tableWithOutput#"NonPrimeVectors" = NonPrimeVectors;
tableWithOutput#"PrimeVectors" = PrimeVectors;
tableWithOutput#"NonBinomialIdeals" = NonBinomialIdeals;
tableWithOutput#"NonPrimeCones" = NonPrimeCones;
tableWithOutput#"PrimeCones" = PrimeCones;
return (#PrimeCones,tableWithOutput)
)


--input list of vectors and starting ideal  of the variety
--output number of good cones and created a mutable hash table with bad vectors, bad ideals, goodvectors, goodideals,verybadvectors,



--computes  all initial ideals 

findInitialIdeals=method(Options => {
	minConvention=>true});

findInitialIdeals(List,Ideal):=opt->(L,I)->(
allIdeals:={};
initialId:=ideal();
i:=0;
    while (i<#L)  do (
    if  opt.minConvention==false then
    (initialId = initialIdeal(L_i, I))
    else (initialId = initialIdeal(-L_i, I));
    allIdeals=append(allIdeals,initialId);
i=i+1;

);
return allIdeals
)

--input the list of weight vectors and the ideal
--output all the initial ideals


--computes the missing binomials of a list of ideals
--input a list of (non-)prime ideals, output a list of missing binomials

findMissingBinomials=method();

findInitialIdeals(Ideal):=(I)->(
j:=0;
k:=0;
PrimaryDec:={};
BinomialComponent:={};
MissingBinomials:={};
PrimaryDec=primaryDecomposition I;
    for j from 0 to #PrimaryDec-1 do(
	if isBinomialIdeal(PrimaryDec)==true then(
	   BinomialComponent=flatten entries(mingens PrimaryDec);
	    for k from 0 to #BinomialComponent-1 do(
        	if isSubset(sub(ideal(BinomialComponent_k),ring I),I)==false then
		MissingBinomials=append(MissingBinomials,BinomialComponent_k);
        
		)
	    )
	else null
	
	)
    ;
return MissingBinomials   
)

findMissingBinomials(List):=(L)->(
i:=0;
j:=0;
k:=0;
PrimaryDec:={};
BinomialComponent:={};
MissingBinomials:={};
for i from 0 to #L-1 do (
    PrimaryDec=primaryDecomposition L_i;
    for j from 0 to #PrimaryDec-1 do(
	if isBinomialIdeal(PrimaryDec_j)==true then(
	   -- print PrimaryDec_j;
	   BinomialComponent=flatten entries(mingens PrimaryDec_j);
	    for k from 0 to #BinomialComponent-1 do(
        
		if isSubset(sub(ideal(BinomialComponent_k),ring L_i),L_i)==false then
		MissingBinomials=append(MissingBinomials,BinomialComponent_k);
        
		)
	    )
	else null
	
	)
    ) ;
return MissingBinomials     
    )

--computes the ideal I' as defined in the paper
--computeModifiedIdeal=(L,V,I)->(
--i:=0;
--R:=QQ[gens ring I,V];

--ListOfRelations:=();
--for i from 0 to #L-1 do(
--	ListOfRelations=append(ListOfRelations,ideal(-V_i)+sub(ideal(L_i),R));
--	);
  --  return (R,sub(I,R)+ideal(ListOfRelations));
   -- )





--computes the map between trop(I) and trop(I'), note that if there is a cone whose projections intersects more maximal cones in their interior then 
--in the output the cone is mapped to nothing

mapMaximalCones=(newR,newM,newL,R4,M4,L4)->(
     --creates the list of the projected rays, forgetting the last two coordinates;
     i:=0;
     j:=0;
     --n is the ambient dimension of the old tropcial variety
     n:=#R4_0;
     mylistOfRays:={};
     while i<#newR do(
     mylistOfRays=append(mylistOfRays,submatrix (matrix {newR_i},{0..n-1}));
     i=i+1
     );
     --creates a list the maximal cones of the original tropical 
     --variety but each cone is represented by the facets
     i=0;
     s:=0;
     c:=0;
     mcone:={};
     mylistOfMaxCones:={};
     while i<#M4 do(
     	 mcone={};
	 s=#M4_0-1;
     	 for j from 0 to s
	  do(c=M4_i_j;mcone=append(mcone,R4_c));  
     	 mylistOfMaxCones=append(mylistOfMaxCones,fourierMotzkin(transpose matrix mcone ,transpose matrix L4));
     	 i=i+1
     	 );
  --looks for a cone of original tropical variety where  newM_i projects
  i=0;
     projCone:={}; 
     projCones:={};
     mcone={};
     
     while i<#newM do(
	 
	 s=newM_i_0;
	 mcone=mylistOfRays_s;
	 for j from 1 to #newM_i-1 do(c=newM_i_j;mcone=mcone||mylistOfRays_c);
	 
     	 j=0;
     	 projCone={"MaxConeNumber",i,newM_i,"=>"};
     	 while(j<#M4)do(
	     
     	     if max(flatten entries (mcone*(mylistOfMaxCones_j_0))   )<=0
      	     and mcone*(mylistOfMaxCones_j_1)==0 
     	     then (
	 	 projCone=append(projCone,M4_j);
		 projCone=append(projCone,"max");
		 projCone=append(projCone,j)
	 )
     ;j=j+1
     );
     if projCone==={"MaxConeNumber",i,newM_i,"=>"} then print i,
     projCones=append(projCones,projCone);
     i=i+1
     );
     return projCones
     )
   
 
 

-------------
--functions to access the hash table computed in findPrimeCones
-------


 

getNonPrimeCones = method(TypicalValue => List)

getNonPrimeCones (MutableHashTable):= M->( M#"NonPrimeCones");

  
getPrimeCones = method(TypicalValue => List)

getPrimeCones (MutableHashTable):= M->( M#"PrimeCones");


getToricIdeals = method(TypicalValue => List)

getToricIdeals (MutableHashTable):= M->( M#"ToricIdeals");


getNonPrimeIdeals = method(TypicalValue => List)

getNonPrimeIdeals(MutableHashTable):= M->( M#"NonPrimeIdeals");


getPrimeVectors = method(TypicalValue => List)

getPrimeVectors (MutableHashTable):= M->( M#"PrimeVectors");

getNonPrimeVectors = method(TypicalValue => List)

getNonPrimeVectors (MutableHashTable):= M->( M#"NonPrimeVectors");


getNonBinomialIdeals = method(TypicalValue => List)

getNonBinomialIdeals (MutableHashTable):= M->( M#"NonBinomialIdeals");












 
 
 
 
-------------------------------------------------------------------------------
-- DOCUMENTATION
------------------------------------------------------------------------------




--"mapMaximalCones",

--"computeModifiedIdeal"
beginDocumentation()
doc ///
    Key
        ToricDegenerations
    Headline
    	a package to compute toric degenerations from the tropicalization of a projective variety
    Description
    	Text
	    This package contains all the functions used 
	     in  "Computing toric degenerations of Flag varieties"[L.Bossinger,S.Lamboglia,K.Mincheva,
		 F.Moammadi]
	    
///
doc ///
    Key 
        getPrimeCones
    Headline
    	
    Usage
        getPrimeCones(H)
    Inputs
	H:HashTable
	    it is one of the outputs of findToricDegenerations	    
    Outputs
    	L:  List
	     of List, each represents a maximal cone 
    Description
    	Text
	  The input is the HashTable which is the output  of  findToricDegenerations.
	  It gives the maximal cones of trop(V(I)) which lead to prime initial ideals.
	  I is  given as input to findToricDegenerations.
///
doc ///
    Key 
        getNonPrimeCones
    Headline
    	
    Usage
        getNonPrimeCones(H)
    Inputs
	H:HashTable
	    it is one of the outputs of findToricDegenerations	    
    Outputs
    	L:  List
	     of List, each represents a maximal cone 
    Description
    	Text
	  The input is the HashTable which is the output  of  findToricDegenerations.
	  It gives the maximal cones of trop(V(I)) which lead to non prime initial ideals.
	  I is  given as input to findToricDegenerations.
///
doc ///
    Key 
        getNonBinomialIdeals
    Headline
    	
    Usage
        getNonBinomialIdeals(H)
    Inputs
	H:HashTable
	    it is one of the outputs of findToricDegenerations	    
    Outputs
    	L:  List
	     of ideals
    Description
    	Text
	  The input is the HashTable which is the output  of  findToricDegenerations.
	  It gives initial ideals corresponding to the maximal cones of trop(V(I)) that are not binomial.
	  I is  given as input to findToricDegenerations.
///
doc ///
    Key 
        getNonPrimeVectors
    Headline
    	
    Usage
        getNonPrimeVectors(H)
    Inputs
	H:HashTable
	    it is one of the outputs of findToricDegenerations	    
    Outputs
    	L:  List
	     of vectors
    Description
    	Text
    	  The input is the HashTable which is the output  of  findToricDegenerations.
    	  It gives the weight vectors w in the maximal cones of trop(V(I)) such that in_w(I) is not prime.
	  I is  given as input to findToricDegenerations.
	    ///

doc ///
    Key 
        getPrimeVectors
    Headline
    	
    Usage
        getPrimeVectors(H)
    Inputs
	H:HashTable
	    it is one of the outputs of findToricDegenerations	    
    Outputs
    	L:  List
	     of vectors
    Description
    	Text
	  The input is the HashTable which is the output  of  findToricDegenerations.
    	  It gives the weight vectors w in the maximal cones of trop(V(I)) such that in_w(I) is prime.
	  I is  given as input to findToricDegenerations.
///
doc ///
    Key 
        getNonPrimeIdeals
    Headline
    	
    Usage
        getNonPrimeIdeals(H)
    Inputs
	H:HashTable
	    it is one of the outputs of findToricDegenerations	    
    Outputs
    	L:  List
	     of ideals
    Description
    	Text
	  The input is the HashTable which is the output  of  findToricDegenerations.
    	  It gives the non-prime initial ideals corresponding to
	  the maximal cones of trop(V(I)). I is given as input to findToricDegenerations.


///
doc ///
    Key 
       getToricIdeals
    Headline
    	
    Usage
        getToricIdeals(H)
    Inputs
	H:HashTable
	    it is one of the outputs of findToricDegenerations	    
    Outputs
    	L:  List
	     of ideals
    Description
    	Text
	  The input is the HashTable which is the output  of  findToricDegenerations.
    	  It gives the toric ideals corresponding to
	  the maximal cones of trop(V(I)). I is given as input to findToricDegenerations.
///

doc ///
    Key 
        findInitialIdeals
    Headline
    	computes the  initial ideals with weight vectors in L
    Usage
        findInitialIdeals(L,I)
    Inputs
	L:List 
	    of weight vectors
	I:Ideal
	    
    Outputs
    	L:  List
	     of initial ideals
    Description
    	Text
    	  It computes the initial ideals with respect to the weight vectors in L. If the list L
	  is the output of  computeWeigthVector(T), where T=trop(V(I)), it computes all 
	   the initial ideals associated to T, hence they are all monomial
	  free. 
///
doc ///
    Key 
        findInitialIdeals
    Headline
    	computes the  initial ideals with weight vectors in L
    Usage
        findInitialIdeals(L,I)
    Inputs
	L:List 
	    of weight vectors
	I:Ideal
	    
    Outputs
    	L:  List
	     of initial ideals
    Description
    	Text
    	  It computes the initial ideals with respect to the weight vectors in L. If the list L
	  is the output of  computeWeigthVector(T), where T=trop(V(I)), it computes all 
	   the initial ideals associated to T, hence they are all monomial
	  free. 
///
doc ///
    Key 
        isBinomialIdeal
    Headline
    	function to check whether an ideal is generated by binomials
    Usage
       isBinomialIdeal(I)
    Inputs
	I:Ideal
    Outputs
    	B: Boolean
	
    Description
    	Text
    	  It checks whether a given ideal is generated by binomials.
///	  
doc ///
    Key 
    	computeWeightVectors
    Headline
    	a function that computes all weight vectors for  initial ideals in the variety
    Usage
        computeWeightVectors(T)
	computeWeightVectors(C,R,L)
    Inputs
        C:List 
	    Maximal Cones of the tropical variety trop(V(I))
	R:List 
	    Rays of the  tropical variety  trop(V(I))
	L:List 
	    Basis of the lineality space of tropical variety  trop(V(I))
	T:
	    a TropicalCycle output of tropicalVariety(I)
	minConvention=> Boolean
	    which indicates the max or min convention.
    Outputs
    	L:  List
	      of weight vectors, one in the relative interior of each maximal cone 
    Description
    	Text
    	  For each cone it computes a weight vectors in the relative interor. By default the all
	  the entries of each vector are made negative, in case of minConvention=>false, then  the entries 
	  will all be positive. This is due to the fact that these vectors will be the input
	  for  initialIdeal.
	  
	  
///
doc ///
    Key 
    	findToricDegenerations
    Headline
    	a function that computes all Groebner toric degenerations of an ideal
    Usage
        findToricDegenerations(L,I)
    Inputs
       
	L:List 
	    of all  weight vectors associated to trop(V(I))
	I:Ideal
	   
	minConvention=> Boolean
	    which indicates the max or min convention.
    Outputs
    	n:List
	    the number of Groebner toric degenerations found
	H:HashTable
	      containing all the data on the Groebner degenerations of I.
    Description
    	Text
    	  It looks for the initial ideals of I that are binomial and prime and outputs the number of 
	  them. It computes also the other non prime initial ideals and checks whether they are binomial
	  or not. It stores in an HashTable all the initial ideals( grouping them in ToricIdeals,NonPrimeIdeals
	  ,nonBinomialIdeals),all the corresponding weight vectors(PrimeVectors,NonPrimeVectors) and
          maximal cones(PrimeCones,NonPrimeCones) of trop(V(I)).
	  The input list L can be computed with computeWeightVectors().
///	  
doc ///
    Key 
    	mapMaximalCones
    Headline
    	
    Usage
        mapMaximalCones(newR,newM,newL,R,M,L)
    Inputs
       
	newR:List 
	       Rays of T'
	newL:List
	       Basis of the lineality space of T'
	newM:List
	       Maximal cones of T'
	R:   List
	       Rays of T
	M:   List
	       Maximal cones  of T
	L:   List
	       Basis of the lineality space of T
    Outputs
    	K:   List
	       of maximal cones of T' with the corresponding image in T
    Description
    	Text 
	  It computes the map \pi:T'\to T  where T=trop (V(I)) for an ideal I and
	  T'=trop(V( I')) where I' is the 
	  ideal described in section 4 of Computing Toric Degenerations of Flag Varieties.
    	  
///
doc ///
    Key 
       findMissingBinomials
    Headline
    	
    Usage
        findMissingBinomials(I)
	findMissingBinomials(L)
    Inputs
       
	J:Ideal
	     
	L:List
	     
    Outputs
        K:List
	       of binomials
    Description
    	Text 
	  The input ideal is assumed to be an initial ideal J=in_w(I)  with respect to a weight vector w that lies in
	  a  maximal cone of trop(V(I)). Moreover the maximal cone has to have multiplicity 1. This
	  implies that in the primary decomposition of J there is only one binomial prime ideal.
	  It computes binomials {f1,..,fs}such that I+(f1,...,fs) is a prime binomial ideal. 
          In case the input is a list of binomial ideals, it does the same for any ideal in the list.
    	  
///

TEST ///
    assert (1+1==2) 
///
