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
"getPrimeIdeals" ,
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
"isBinomialIdeal",
"computeModifiedIdeal"
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


--computes the missing binomials
--input a list of (non-)prime ideals, output a list of missing binomials
findMissingBinomials=(L)->(
i:=0;
j:=0;
k:=0;
PrimaryDec:={};
MissingBinomials:={};
for i from 0 to #L-1 do (
    PrimaryDec=primaryDecomposition L_i;
    for j from 0 to #PrimaryDec-1 do(
	if isBinomialIdeal(PrimaryDec_j)==true then(
	   -- print PrimaryDec_j;
	    for k from 0 to #entries(mingens PrimaryDec_j)-1 do(
		--print (entries(mingens PrimaryDec_j))_k;
		if isSubset(ideal((entries(mingens PrimaryDec_j))_k),L_i)==false then
		MissingBinomials=append(MissingBinomials,(entries(mingens PrimaryDec_j))_k);
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

mapMaximalCones=(newR,newM,newL,R,M,L)->(
     --creates the list of the projected rays, forgetting the last two coordinates;
     i:=0;
     j:=0;
     --n is the ambient dimension of the old tropcial variety
     n:=#R_0;
     mylistOfRays:={};
     while i<#newR do(
     mylistOfRays=append(mylistOfRays,submatrix (matrix {newR_i},{0..n-1}));
     i=i+1
     );
     print "ok";
   
     --creates a list the maximal cones of the original tropical 
     --variety but each cone is represented by the facets
     i=0;
     mcone:={};
     mylistOfMaxCones:={};
     while i<#M do(
     	 provList:={};
     	 mcone={};
     	 for j from 0 to #M_0-1 do(mcone=append{R_M_j});     
     	 mylistOfMaxCones=append(mylistOfMaxCones,fourierMotzkin(transpose matrix mcone ,transpose matrix L));
     	 i=i+1
     	 );
     print "done";
  --looks for a cone of original tropical variety where  newM_i projects
     i=0;
     projCone:={}; 
     projCones:={};
     mcone={};
     while i<#newM do(
	 mcone=mylistOfRays_newM_0;
	 for j from 1 to #newM-1 do(mcone=mcone||mylistOfRays_newM_j);
     	 j=0;
     	 projCone={"MaxConeNumber",i,newM_j,"=>"};
     	 while(j<#M)do(
     	     if max(flatten entries (mcone*(mylistOfMaxCones_j_0))   )<=0
      	     and mcone*(mylistOfMaxCones_j_1)==0 
     	     then (
	 	 projCone=append(projCone,M_j);
		 projCone=append(projCone,"max");
		 projCone=append(projCone,j)
	 )
     ;j=j+1
     );
     if projCone==={"MaxConeNumber",i,newM_j,"=>"} then print i;
     projCones=append(projCones,projCone);
     i=i+1
     );
     return projCones
     );
 
 
 

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
beginDocumentation()
doc ///
    Key
        ToricDegenerations
    Headline
    	a package to compute toric degenerations from the tropicalization of a projective variety
    Description
    	Text
	    This package contains all the functions used  in  "Computing toric degenerations of Flag varieties"[L.Boss]
	    
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
	R:List
	L:List
	T:Table
    Outputs
    	L:  List
    Description
    	Text
    	  For each cone it computes a weight vectors in the relative interor
	  
///	  
TEST ///
    assert (1+1==2) 
///
