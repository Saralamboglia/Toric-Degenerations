--constructing the ideal of Flag 4

R=QQ[p_(0),p_(1),p_(2),p_(3),p_(0,1),p_(0,2),p_(0,3),p_(1,2),p_(1,3),p_(2,3),p_(0,1,2),p_(0,1,3),p_(0,2,3),p_(1,2,3),a,b,c,d,e,f,g,h]

M=matrix{{p_(0),p_(1),p_(2),p_(3)},{a,b,c,d}};

minorsM=minors(2,M);

I1=ideal();
indexs=subsets(4,2);for i from 0 to 5 do (oneindex=toSequence indexs_i;I1=I1+ideal(p_oneindex-(minorsM_i));i=i+1);I1

M=matrix{{p_(0),p_(1),p_(2),p_(3)},{a,b,c,d},{e,f,g,h}};

minorsM=minors(3,M);

I2=ideal();
indexs=subsets(4,3);for i from 0 to 3 do (oneindex=toSequence indexs_i;I2=I2+ideal(p_oneindex-(minorsM_i));i=i+1);I2

Iflag4= eliminate({a,b,c,d,e,f,g,h},I1+I2);

-- tropicalizing it

R2=QQ[p_(0),p_(1),p_(2),p_(3),p_(0,1),p_(0,2),p_(0,3),p_(1,2),p_(1,3),p_(2,3),p_(0,1,2),p_(0,1,3),p_(0,2,3),p_(1,2,3)]

Iflag4=sub(Iflag4,R2);

Tflag4=tropicalVariety(Iflag4,true);

raysTflag4=getRays Tflag4;

maximalConesTflag4=getMaximalCones Tflag4;

multiplicitiesTflag4=getMultiplicities Tflag4;

linealitySpaceTflag4=getLinealitySpace Tflag4;




--constructing the ideal of Flag 5

R1=QQ[a,b,c,d,e,x,y,z,w,u,k,q,r,s,t,p_0..p_4,apply(subsets(5,2),i->p_i),apply(subsets(5,3),i->p_i),apply(subsets(5,4),i->p_i)]

M=matrix{{p_0,p_1,p_2,p_3,p_4},{a,b,c,d,e}};

minorsM=minors(2,M);

I1=ideal();
indexs=subsets(5,2);for i from 0 to 9 do (I1=I1+ideal(sub(p_(indexs_i),R1)-(minorsM_i));i=i+1);I1

M=matrix{{p_0,p_1,p_2,p_3,p_4},{a,b,c,d,e},{x,y,z,w,u}};

minorsM=minors(3,M);

I2=ideal();
indexs=subsets(5,3);for i from 0 to 9 do (I2=I2+ideal(sub(p_(indexs_i),R1)-(minorsM_i));i=i+1);I2


M=matrix{{p_0,p_1,p_2,p_3,p_4},{a,b,c,d,e},{x,y,z,w,u},{k,q,r,s,t}};

minorsM=minors(4,M);

I3=ideal();
indexs=subsets(5,4);for i from 0 to 4 do (I3=I3+ideal(sub(p_(indexs_i),R1)-(minorsM_i));i=i+1);I3

Iflag5=eliminate({a,b,c,d,e,x,y,z,w,u,k,q,r,s,t},I1+I2+I3);

isPrime Iflag5;--note that takes too long 













 











