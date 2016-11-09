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
export {"computeModifiedIdeal",
"computeWeightVectors",
"getNonPrimeIdeals",
"getToricIdeals",
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
"isBinomialIdeal",
"scanInitialIdeals"}

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
scanInitialIdeals=method(Options => {
	minConvention=>true});






scanInitialIdeals(List,Ideal):=opt->(L,I)->(
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
	    print "Found One Prime";
	    )
     	else ( 
	    NonPrimeIdeals=append(NonPrimeIdeals,initialId);
            NonPrimeVectors=append(NonPrimeVectors,L_i);
            NonPrimeCones=append(NonPrimeCones,i);
	    print "Found One Non Prime";
	    
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

scanInitialIdeals(Ideal):=opt->(I)->(
    if opt.minConvention==false then print "Tropicalize it with gfan and then use
   scanInitialIdeals(computeWeightVectors(T,minConvention=>false),I,minConvention=>false)"
   else(
       T:=tropicalVariety(I);
    scanInitialIdeals(computeWeightVectors(T),I);
    )
    )

scanInitialIdeals(TropicalCycle,Ideal):=opt->(T,I)->(
    scanInitialIdeals(computeWeightVectors(T),I);
    )
    

--input list of vectors and starting ideal  of the variety
--output number of good cones and created a mutable hash table with bad vectors, bad ideals, goodvectors, goodideals,verybadvectors,



--computes  all initial ideals 

findInitialIdeals=method(Options => {
	minConvention=>true});
findInitialIdeals(TropicalCycle,Ideal):=opt->(T,I)->(
    findInitialIdeals(computeWeightVectors(T,opt.minConvention),I)
    )
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

findMissingBinomials(Ideal):=(I)->(
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
--L is the list of binomials,V is the set of new variables ,H list of homogenizing variables,I old ideal
computeModifiedIdeal=(L,V,H,I)->(
i:=0;
R:=QQ[gens ring I,V,H];
n:=#(gens ring I);
k:=#V;
S:=0;
A:=gens R;
J:=sub(I,R);

for i from 0 to #L-1 do(
	J=homogenize(J+ideal((gens R)_{n+i}-{sub(L_i,R)}),A_(n+k+i));
	);
  return sub(J,R);
  )
 


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

--functions that checks if there are toric Groebner degenerations, if there are it gives them,
--otherwise it re-embedds the variety and checks if from the new embedding it is
--possible to find new toric degenerations(equivalently if there are prime cones)
--the ideal I has not to contain monomials
findToricDegenerations=method(Options => {
	minConvention=>true});


findToricDegenerations(Ideal):=opt->(I)->(
   if opt.minConvention==false then print "Tropicalize it with gfan and then use
   ToricDegenerations(T)"
   else(
       if isHomogeneous(I)==false then print "The ideal I is not homogeneous,please homegenize before
       using this function"
       else(
       	   T:=tropicalVariety(I);
	  
	   ToricOutput:=scanInitialIdeals(computeWeightVectors(T),I);
	   if (ToricOutput)_0>0
	   then return ToricOutput 
       	   else ( 
	       print "There are no toric degenerations, Re-embedding in process";
	       MissingBinomials:=findMissingBinomials(getNonPrimeIdeals(ToricOutput_1));
	       print MissingBinomials;
	       n:=#MissingBinomials;
	       Y:=symbol Y;
	       H:=symbol H;
	       I':=computeModifiedIdeal(MissingBinomials,toList(Y_1..Y_n),toList(H_1..H_n),I);
	       print toString I';
	       
	       --there is a trick since gfan starting cone does not work,hence I use 
	       --gfan_bruteforce pretending that I' is not prime
	 
	       T':=tropicalVariety(I',Prime=>false);
	       ToricOutput':=scanInitialIdeals(computeWeightVectors(T'),I');
	       return (ToricOutput',T')
	           
    )
    )
    )
    )

 

-------------
--functions to access the hash table computed in findPrimeCones
-------


 

getNonPrimeCones = method();


getNonPrimeCones (MutableHashTable):= M->( M#"NonPrimeCones");

  
getPrimeCones = method();


getPrimeCones (MutableHashTable):= M->( M#"PrimeCones");


getToricIdeals = method();


getToricIdeals (MutableHashTable):= M->( M#"ToricIdeals");


getNonPrimeIdeals = method();


getNonPrimeIdeals(MutableHashTable):= M->( M#"NonPrimeIdeals");


getPrimeVectors = method();


getPrimeVectors (MutableHashTable):= M->( M#"PrimeVectors");

getNonPrimeVectors = method();


getNonPrimeVectors (MutableHashTable):= M->( M#"NonPrimeVectors");


getNonBinomialIdeals = method();


getNonBinomialIdeals (MutableHashTable):= M->( M#"NonBinomialIdeals");












 
 
 
 
-------------------------------------------------------------------------------
-- DOCUMENTATION
------------------------------------------------------------------------------



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
		 F.Moammadi]. 
	    
///
doc ///
    Key 
        getPrimeCones
    Headline
    	
    Usage
        getPrimeCones(H)
    Inputs
	H:HashTable
	    it is one of the outputs of  scanInitialIdeals	    
    Outputs
    	L:  List
	     of List, each represents a maximal cone 
    Description
    	Text
	  The input is the HashTable which is the output  of  findToricDegenerations.
	  It gives the maximal cones of trop(V(I)) which lead to prime initial ideals.
	  I is  given as input to scanInitialIdeals.
///
doc ///
    Key 
        getNonPrimeCones
    Headline
    	
    Usage
        getNonPrimeCones(H)
    Inputs
	H:HashTable
	    it is one of the outputs of  scanInitialIdeals    
    Outputs
    	L:  List
	     of List, each represents a maximal cone 
    Description
    	Text
	  The input is the HashTable which is the output  of  findToricDegenerations.
	  It gives the maximal cones of trop(V(I)) which lead to non prime initial ideals.
	  I is  given as input to  scanInitialIdeals.
///
doc ///
    Key 
        getNonBinomialIdeals
    Headline
    	
    Usage
        getNonBinomialIdeals(H)
    Inputs
	H:HashTable
	    it is one of the outputs of  scanInitialIdeals	    
    Outputs
    	L:  List
	     of ideals
    Description
    	Text
	  The input is the HashTable which is the output  of  findToricDegenerations.
	  It gives initial ideals corresponding to the maximal cones of trop(V(I)) that are not binomial.
	  I is  given as input to  scanInitialIdeals.
///
doc ///
    Key 
        getNonPrimeVectors
    Headline
    	
    Usage
        getNonPrimeVectors(H)
    Inputs
	H:HashTable
	    it is one of the outputs of  scanInitialIdeals	    
    Outputs
    	L:  List
	     of vectors
    Description
    	Text
    	  The input is the HashTable which is the output  of  findToricDegenerations.
    	  It gives the weight vectors w in the maximal cones of trop(V(I)) such that in_w(I) is not prime.
	  I is  given as input to  scanInitialIdeals.
///

doc ///
    Key 
        getPrimeVectors
    Headline
    	
    Usage
        getPrimeVectors(H)
    Inputs
	H:HashTable
	    it is one of the outputs of  scanInitialIdeals 
    Outputs
    	L:  List
	     of vectors
    Description
    	Text
	  The input is the HashTable which is the output  of  findToricDegenerations.
    	  It gives the weight vectors w in the maximal cones of trop(V(I)) such that in_w(I) is prime.
	  I is  given as input to  scanInitialIdeals.
///
doc ///
    Key 
        getNonPrimeIdeals
    Headline
    	
    Usage
        getNonPrimeIdeals(H)
    Inputs
	H:HashTable
	    it is one of the outputs of  scanInitialIdeals    
    Outputs
    	L:  List
	     of ideals
    Description
    	Text
	  The input is the HashTable which is the output  of  findToricDegenerations.
    	  It gives the non-prime initial ideals corresponding to
	  the maximal cones of trop(V(I)). I is given as input to  scanInitialIdeals.


///

doc ///
    Key 
       getToricIdeals
    Headline
    	
    Usage
        getToricIdeals(H)
    Inputs
	H:HashTable
	    it is one of the outputs of  scanInitialIdeals	    
    Outputs
    	L:  List
	     of ideals
    Description
    	Text
	  The input is the HashTable which is the output  of   scanInitialIdeals.
    	  It gives the toric ideals corresponding to
	  the maximal cones of trop(V(I)). I is given as input to  scanInitialIdeals.
///

doc ///
    Key 
        findInitialIdeals
    Headline
    	computes the  initial ideals with weight vectors in L
    Usage
        findInitialIdeals(L,I)
	findInitialIdeals(T,I)
    Inputs
	L:List 
	    of weight vectors
	I:Ideal
	T:
	   TropicalCycle, which is trop()
	    
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
    	scanInitialIdeals
    Headline
    	a function that computes all Groebner toric degenerations of an ideal
    Usage
        scanInitialIdeals(L,I)
	scanInitialIdeals(T,I)
	scanInitialIdeals(I)
    Inputs
       
	L:List 
	    of all  weight vectors associated to trop(V(I))
	I:Ideal
	T:
	   TropicalCycle, which is trop(V(I))
	   
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
        findToricDegenerations
    Headline
    	a function that computes   toric degenerations of a variety V(I)
    Usage
        findToricDegenerations(I)
    Inputs
       
	I:Ideal
	   
	minConvention=> Boolean
	    which indicates the max or min convention.
    Outputs
    	n:List
	    the number of Groebner toric degenerations found
	H:HashTable
	      containing all the data on the Groebner degenerations of I ot I'
	T:
	  TropicalCycle which is  trop(V(I)) or trop(V(I'))
    Description
    	Text
    	  It first applies scanInitialIdeals to I to check whether there  prime(hence binomial)
	   initial ideals. If there is at least one is one it outputs the number of them,
	  the HashTable created by scannInitialIdeals and the tropicalization trop(V(I)).
	  If there are no prime initial ideals then it considers I' computed by computeModifiedIdeal
	  (see section 4 of the paper). It applies again scanInitialIdeals to see if there are 
	  prime initial ideals of I'. In this case the output will be the number of prime initial ideals
	  of I', the HashTable created by scannInitialIdeals and the tropicalization trop(V(I')).
	  Note that as in the case of computeModifiedIdeal, the ideal I' has new variables of the form
	  Y_i, H_i, hence these should not be used in the input ideal.
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
	  ideal described in section 4 of the paper.
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
	    of Ideals 
    Outputs
        K:List
	       of binomials
    Description
    	Text 
	  The input ideal is assumed to be an initial ideal J=in_w(I)  with respect to a weight vector w that lies in
	  a  maximal cone of trop(V(I)). Moreover the maximal cone has to have multiplicity 1. This
	  implies that in the primary decomposition of J there is only one binomial prime ideal.
	  It computes binomials {f1,..,fs}such that J+(f1,...,fs) is a prime binomial ideal. 
          In case the input is a list of binomial ideals, it does the same for any ideal in the list.
    	  
///


doc ///
    Key 
       computeModifiedIdeal
    Headline
    	
    Usage
        computeModifiedIdeal(L,V,H,I)
    Inputs
	I:Ideal
	     
	L:List
	    of binomials to add
        V:List
	    of new variables
        H:List
	   of new variables	     
    Outputs
        I':Ideal
	       
    Description
    	Text 
	  Note that the ideal I' has new variables of the form
	  Y_i, H_i, hence these should not be used in the input ideal.
	  The  ideal I is assumed to be  prime and homogeneous. It first computes I'=I+(V_0-L_0) and it 
	  homogenizes it with respect to H_0. Then it does the same starting with I',V_1,L_1,H_1 and analogously
	  for the rest of L_i,V_i. If I is an ideal in QQ[x_1,...,x_n] then  the ideal I'
	  considered in each step is contained in
	  a new ring with variables giving x_1,...,x_n and by the elements in the lists V and H.
	  The output is the ideal I' described in section 4 of the paper. Note that this is a homogeneous and prime  ideal(if I is prime).
	  The list of binomials can be computed with findMissingBinomials.
	  	 
    	  
///



TEST ///
    assert (1+1==2) 
///
