(* ::Package:: *)

BeginPackage["NMStools`"]


(*Main functions that will be used for eigenvalues and non-modal calculations*)
ClearAll[EigenNDSolve,pseudoModes,plotEValues,checkSolutions,shootNDsolve,generalChebyMat,
	constructDmat,CCapprox, transCo,eigsToInterp,funcFromCheby,myTiming,intFuncsPlot]

EigenNDSolve::usage="EigenNDsolve[eqns, vars, range, evalvar, Options]
Solves the list of eqns for the eigenfunctions and eigenvalues using a spectral expansion
in Chebyshev polynomials.
vars is a list of the dependent variables, 
range is a list with {indepentent variable, min, max} (as for NDSolve),
evalvar is the eigenvalue parameter (must be linear in this).

Options:
Npoly\[Rule]80 - Number of polynomials to use in Chebyshev expansion (actually Npoly+1 due to \!\(\*SubscriptBox[\(T\), \(0\)]\))
returnInterp\[Rule]True - If True, returns the eigenfunctions as InterpolatingFunction objects.
	If False, returns a list of the Chebyshev coefficients for each eigenstate.
noRowReduce\[Rule]False - If True, over-rides solving for the boundary conditions and forming a
	reduced matrix.
precision\[Rule]MachinePrecision - Precision of the calculation. If not MachinePrecision,
	calculation will be much slower; however this is necessary for accurate calculation
	with certain types of problems.
returnMatrix->False - With normal inputs, returns a the Chebyshev matrix. If the output is 
	then used as eqns, it skips calculating the matrix. Useful for parameter scans."

pseudoModes::usage="pseudoModes[esys, K, inP, Options]
Returns functions to calculate the spatial structure and growth of the fastest growing 
pseudo-modes at a given time.
esys is the list of eigenvalues and eigenvectors.
K is the number of eigenvectors to use in calculating fastest growing mode.
inP is the inner product to use. This can be either a matrix, in which case the inner
	product is assumed to be v\[ConjugateTranspose].M.v, or a function taking arguments of the form 
	inP[{{\!\(\*SubscriptBox[\(\[Xi]\), \(1\)]\),\!\(\*SubscriptBox[\(\[Xi]\), \(2\)]\)},{\!\(\*SubscriptBox[\(\[Eta]\), \(1\)]\),\!\(\*SubscriptBox[\(\[Eta]\), \(2\)]\)},...}] where {\[Xi], \[Eta],...} are the dependent variables.
Output functions, {G[t],\[Xi][t]}
G[t] is the maximum growth at a given time
\[Xi][t,\!\(\*SubscriptBox[\(t\), \(0\)]\)] is the perturbation at time t that maximizes the growth at \!\(\*SubscriptBox[\(t\), \(0\)]\)

Options:
cholPrecision\[Rule]MachinePrecision - Precision used in calculating the Cholesky factorizations
	and Hermitian matrices. May need to be increased if PositiveDefinite warning encountered.
HermitianFix\[Epsilon]\[Rule]\!\(\*SuperscriptBox[\(10\), \(-14\)]\) - Adds \[Epsilon] DiagonalMatrix to matrix before taking CholeskyFactorization, to 
	help make it positive definite. Safer to increase cholPrecision. 
removeLarEvalsCutoff\[Rule]\!\(\*SuperscriptBox[\(10\), \(3\)]\) - Automatically remove any eigenvalues that are larger than this."

plotEValues::usage="plotEValues[evalues, range, Options]
Plots eigenvalues in a convenient way, displaying their value and number with a mouse-over.
range is a nested list {\!\(\*SubscriptBox[\(x\), \(range\)]\), \!\(\*SubscriptBox[\(y\), \(range\)]\)}.

Options:
estateList->False - List of eigen-functions to plot when a given point is clicked.";


checkSolutions::usage="checkSolutions[eqns, vars, range, evalvar, solutions]
Similar syntax to EigenNDSolve but range includes y range, evalvar second argument 
is number of the eigenfunction, and solutions is the solution as returned by 
EigenNDSolve.";

shootNDsolve::usage="shootNDsolve[eqns, vars, range, esys ,Options]
Solves eqns using standard NDSolve method, vars and range as for EigenNDSolve.
esys should be in form {evar, number, esystem} where esystem is the output 
of EigenNDSolve (with InterpolatingFunctions).

Options:
useRightBC\[Rule]False - Take boundary conditions for right rather than left boundary. 
evalOffset\[Rule]0. - Offset the eigenvalue by given amount";

constructDmat::usage = "constructDmat[ord_, size_List ]";

generalChebyMat::usage = "generalChebyMat[size_List, zfuncs_List]
Prodcues a Chebyshev matrix of size size, where zfuncs is a list of multiplying each derivative (highest
derivative first).

Options:
precision\[Rule]MachinePrecision";

CCapprox::usage = "CCaprox[func_, n_, eigv_, prec_:MachinePrecision]
Approximate function by Chebyshev sum.
eigv is a symbolic parameter (doesn't matter if there isn't one).";

transCo::usage = "transCo[eqn_, varQ_, xrange_, zrange_]
Trasform an equation to a different variable range. 
varQ is a function that returns True for the dependent variables.";

funcFromCheby::usage ="funcFromCheby[chebco_,rang_:{-1,1}, Options]
Produces interpolating function from list of Chebyshev coefficients, chebco.
Options:
Npoints\[Rule]101 - Number of points over which to interpolate.
returnList\[Rule]False - If True, returns list rather than InterpolationFunction";

eigsToInterp::usage = "eigsToInterp[eigsys_,range_:{-1,1}]
Turns the eigenfunctions in eigsys (as returned by EigenNDSolve) into InterpolatingFunctions.";

myTiming::usage = "myTiming[expr]
Prints out Timing without changing return value";

intFuncsPlot::usage = "intFuncsPlot[intfuncs_List]
Plots list of InterpolatingFunctions in a nice way showing real and imaginary parts";




(*Functions to use for deriving equations*)
ClearAll[linearize]
linearize::usage = "linearize[expr_, vars_List, indvars_List]
Linearizes expr in the variables vars. indvars are the independent variables.";



(*Old functions for use with the finite difference solver*)
ClearAll[FDEigenNDSolve,linearNPDESolve,polyEigen]

FDEigenNDSolve::usage="Similar to EigenNDSolve syntax but uses finite differences. 
Doesn't work well at moderate Reynolds numbers";

polyEigen::usage = "Solves generalized polynomial eigenvalue problem";

linearPDESolve::usage ="Solves the linear PDE in time, doesn't work well at all due to numerical
instability";



(* ::Section:: *)
(*Main package*)


Begin["`Private`"]


(* ::Subsection:: *)
(*Auxilary functions*)


(* ::Subsubsection:: *)
(*Public*)


constructDmat[ord_,size_List]:=Block[{lsiz,arrrules,tempmat},
If[ord==0,
SparseArray[Band[{1,1}]->ConstantArray[1,Min[size]],size],
lsiz=Max[size];
(*Since Dmat will be very common, construct using formula for speed*)
arrrules=Flatten[Table[{{1,2j}->2j-1,{i+1,i+2j}->2(i+2j-1)},{i,1,lsiz},{j,1,lsiz}]]/.({_,n_/;n>lsiz}->_):>Sequence[];
tempmat=SparseArray[arrrules,{lsiz,lsiz}];
MatrixPower[tempmat,ord][[1;;size[[1]],1;;size[[2]]]]]
]



CCapprox[func_,n_,eigv_,prec_:MachinePrecision]:=Block[{cnodes,funclist,cc},
cnodes=N[Cos[Pi Range[0,n]/n],prec];
funclist={#/.eigv->0,Coefficient[#,eigv]}&[func/@cnodes];
cc=FourierDCT[#,1]Sqrt[2/n]&/@funclist;
cc[[All,{1,-1}]]/=2;
cc=SetPrecision[cc,prec];
cc[[1]]+eigv cc[[2]]]


transCo[eqn_,varQ_,xrange_,zrange_]:=
	Module[{x=xrange[[1]],xr=xrange[[2;;]],z=zrange[[1]],zr=zrange[[2;;]],trans},
		trans={Subtract@@xr/Subtract@@zr,(xr[[2]]zr[[1]]-xr[[1]]zr[[2]])/(zr[[1]]-zr[[2]])+Subtract@@xr/Subtract@@zr #&,
				(zr[[2]]xr[[1]]-zr[[1]]xr[[2]])/(xr[[1]]-xr[[2]])+Subtract@@zr/Subtract@@xr #&};
		eqn/.Derivative[n_][f_?varQ][v_]:>1/ (trans[[1]]^n) Derivative[n][f][trans[[3]][v]]/.
			{f_?varQ[v_]:>f[trans[[3]][v]]}/.x->trans[[2]][z]//ExpandAll
]



ClearAll[calc$chebyshevLists]
(*Produces interpolating function from a list of Chebyshev coefficients*)
calc$chebyshevLists[npol_,npoints_]:=$chebyshevLists=Table[N[ChebyshevT[i,Range[-1,1,2/(npoints-1)]]],{i,0,npol}];
calc$chebyshevLists[120,101];
Options[funcFromCheby]={Npoints->101,returnList->False};
funcFromCheby[chebco_,rang_:{-1,1},OptionsPattern[]]:=Module[{},
	If[Length[chebco]>Length@$chebyshevLists||OptionValue[Npoints]!=Length@First@$chebyshevLists,
	calc$chebyshevLists[Length[chebco],OptionValue[Npoints]]];
	If[OptionValue[returnList],
		chebco.$chebyshevLists[[1;;Length[chebco]]],
		ListInterpolation[chebco.$chebyshevLists[[1;;Length[chebco]]],{rang}]
	]
]

eigsToInterp[eigsys_,range_:{-1,1}]:=Block[{esys=eigsys},
esys[[2;;]]=Map[funcFromCheby[#,range]&,esys[[2;;]],{2}];
esys]



Options[shootNDsolve]={useRightBC->False,evalOffset->0.0};
(*Solve equations using NDSolve as a check.*)
(*esys should be in form {evar, number, esystem}, where esystem is output from EigenNDSolve including solution of differential equations as InterpolatingFunctions (it takes the initial conidtions from this)*)
shootNDsolve::bounds="Boundaries not found correctly!";
shootNDsolve[eqns_List,vars_List,range_List,esys_List,OptionsPattern[]]:=
Module[{x,dord,intfuncs,bounds,initconds,lr,
soln,eval},
x=range[[1]];
dord=Max[Cases[eqns,Derivative[n_][#][_]:>n,Infinity]]&/@vars;
(*Find initial conditions and boundaries*)
intfuncs=Rest[esys[[3]]][[All,esys[[2]]]];
lr=If[TrueQ@OptionValue[useRightBC],range[[3]],range[[2]]];
initconds=Flatten[
Function[{n},(D[vars[[n]][x]==intfuncs[[n]][x],{x,#}]&/@Range[0,dord[[n]]-1])/.x->lr]/@Range[Length[vars]]
];

(*Solve differential equation*)
eval=esys[[3,1,esys[[2]]]]+OptionValue[evalOffset];
soln=NDSolve[Join[eqns,initconds]/.esys[[1]]->esys[[3,1,esys[[2]]]],vars,range][[1,All,2]];

Show[intFuncsPlot[soln],PlotLabel->"\[Omega] = "<>ToString[eval]]


]


(*Plots a list of interpolating functions in a nicish way, showing real and imaginary parts*)
intFuncsPlot[intfuncs_List]:=Module[{plotlist,range,collist,x},
plotlist=Flatten[{Re@#,Im@#}&/@Through[intfuncs[x]]];
range=intfuncs[[1,1,1]];

collist=Flatten[{#,Directive[#,Dashed]}&/@(ColorData[1,#]&/@Range[Length[intfuncs]])];

Plot[plotlist,Evaluate[{x,Sequence@@range}],PlotRange->All,PlotStyle->collist]
]




(*Better version of linearize, works on expressions in denominator and multiple variables*)
linearize[expr_,vars_List,indvars_List]:=Module[{dlist,fullvars,t},
(*create full list*)
dlist=Permutations[Flatten[ConstantArray[indvars,4]],4];
fullvars=Flatten[
D[Through[vars@@indvars],Sequence@@#]&/@dlist
];
Expand[Normal[Series[expr/.(#->t #&/@fullvars),{t,0,1}]]/.t->1]
]



SetAttributes[myTiming,HoldAll];
(*Setup to be able to add together multiple loops, could add this from VEST*)
myTiming[f_]:=Module[{outlist},
	$myTimingIter=0;
	outlist=AbsoluteTiming[Reap[f,timingtag]];
	Print[Head[Unevaluated[f]],": ",outlist[[1]] "s"];
	If[Length[outlist[[2,2]]]!=0,
		Print[Flatten[outlist[[2,2]]][[1]]," total: ",Total[Drop[Flatten[outlist[[2,2]]],1]]]];
	outlist[[2,1]]
];


(*Plug solutions given by EigenNDSolve back into equations to check errors*)
(*Simialr syntax to EigenNDSolve but range includes y range, evalvar second argument is number of the eigenfunction, and solutions is the solution as returned by EIgenNDSOlve*)
checkSolutions[eqns_List,vars_List,range_List,evalvar_List,solutions_List]:=
Plot[Evaluate[Flatten[{Re@#,Im@#}&/@
(Subtract@@@(eqns/.Thread[Join[{evalvar[[1]]},vars]->solutions[[All,evalvar[[2]]]]]))]],
Evaluate[{range[[1]],Sequence@@range[[2,1]]}],
PlotRange->range[[2]],PlotStyle->Flatten[{ColorData[1,#],Directive[Dashed,ColorData[1,#]]}&/@Range[Length[eqns]]]
]





Options[generalChebyMat]={precision->MachinePrecision};
generalChebyMat[size_List,zfuncs_List,OptionsPattern[]]:=Module[{T,dt,phivec,chebexpans,diffeq,dum},
	defineAbCheb[T,dt];
	phivec=(T[#]&/@Range[0,size[[2]]]);
	chebexpans=Map[Dot[#,phivec]&,Chop[CCapprox[#,size[[2]],dum,OptionValue[precision]],10^-14]&/@zfuncs];
	diffeq=Total@ParallelMap[ExpandAll, 
		chebexpans  (Nest[ParallelMap[Expand[dt[#]]&,#,$dcon]&,phivec,#]&/@Range[Length[zfuncs]-1,0,-1]),$dcon];
	Join[{diffeq/.T[_]->0},Coefficient[diffeq,T[#]]&/@Range[size[[1]]]]
]


(* ::Subsubsection:: *)
(*Private*)


ClearAll[RowR]
RowR[mat_,{n_,m_}]:=Block[{matt,i},
matt=Table[-(mat[[i,m]]/mat[[n,m]])mat[[n]]+mat[[i]],Evaluate[{i,Complement[Range[Length[mat]],{n}]}]];
Chop[Insert[matt,mat[[n]],n],10^-14]]


ClearAll[checkBCform]
checkBCform[BCs_,vars_,dord_]:=
Block[{totalBCs=Length[BCs],BClist,BCmatch,fornone,RRQ,nBCs,lens,vl,vcl,overlaps,numvars,numBCfin},
BClist=Function[{var},Map[Cases[#,var[_]|Derivative[_][var][_],{0,Infinity}]&,BCs]]/@vars;
If[Or@@Thread[Abs[Flatten[BClist][[All,1]]]!=1],
	Message[EigenNDSolve::bounds,Thread[BCs==0]];Abort[]];
(*Don't want two of the same variable in each BC list, then double count*)
fornone=If[Length[#]!=0,{First@#},#]&;
BClist=Map[fornone,BClist,{2}];
(*Check to see if BCs can be "row reduced" from the matrix*)
nBCs=Length/@(Flatten/@BClist);
BCmatch=Thread[dord<=nBCs];

(*The True/False returned refers to whether the row reduction step should be performed*)
RRQ=If[Length[BClist[[1]]]==Total[dord] && And@@BCmatch,True,
	Message[EigenNDSolve::noVBC];False];

(*Choose the number of extra polynomials needed for each variable*)
lens=Map[Length,BClist,{2}];
vl=Range[Length[vars]];
(*Finds the number of overlapping boundary conditions for each variable*)
overlaps=Function[{vnum},Total[Flatten[If[lens[[vnum,#]]!=0,
	Boole[#!=0]&/@lens[[Complement[vl,{vnum}],#]],0]&/@Range[totalBCs]]]]/@vl;

If[Length[BClist[[1]]]!=Total[dord],
(*If number of boundary conditions doesn't match dord, just reduce some choice of those that overlap*)
	numvars[nbc_,maxred_]:=Block[{i=nbc,j=maxred},While[j>0,--i;--j];i];
	numBCfin=MapThread[numvars,{nBCs,overlaps}],
(*Else work out possible reduction from dord*)
(*number of BCs,possible reduction,dord*)
	numvars[nbc_,maxred_,do_]:=Block[{i=nbc,j=maxred},While[i>do&&j>0,--i;--j];i];
	numBCfin=MapThread[numvars,{nBCs,overlaps,dord}]
];
{numBCfin,RRQ}

]



ClearAll[defineAbCheb]
defineAbCheb[T_,dt_]:=Module[{},
(*Call this new each time to stop variable overlap*)
ClearAll[dt,T];
(*derivative properties (needed if we have a(x)\!\(
\*SubscriptBox[\(\[PartialD]\), \(x\)]\[Xi]\) since only \!\(
\*SubscriptBox[\(\[PartialD]\), \({x, n}\)]\[Xi]\) is stored*)
dt[a_+b_]:=dt[a]+dt[b];
dt[a_ b_]:=a dt[b]+b dt[a];
dt[a_?NumericQ]=0;
(*Derivatives of polynomials*)
T[0]:=1;
dt[T[n_]]:=If[EvenQ[n-1],Expand[n (2Sum[T[j],{j,0,n-1,2}]-1)],Expand[n(2Sum[T[j],{j,1,n-1,2}])]];
(*Products of polynomials*)
T/:T[i_Integer]T[j_Integer]:=1/2(T[i+j]+T[Abs[i-j]]);
T/:T[j_Integer]^2:=1/2(T[2j]+T[0]);
(*Rules for boundaries*)
{T[j_]:>(-1)^j,T[j_]:>1}
]





ClearAll[remEnds,addEnds]
remEnds=Most@Rest@#&;
addEnds=Join[{0},#,{0}]&;



ClearAll[removeLargeEvalues];
(*Remove modes with very large eigenvalues from the list*)
(*Works with either interpolating functions of lists as returned by EigenNDSolve*)
removeLargeEvalues[esys_List,cutoff_:10^5]:=Module[{},
(*Remove grid at end of list*)
esys[[All,First[Position[Thread[Abs[esys[[1]]]<cutoff],True]\[Transpose]]]]
]

ClearAll[removeEvalueList];
(*Remove list of numbers from evalue list*)
removeEvalueList[esys_List,remove_List]:=Module[{len,nevals},
len=Length[esys]-If[Head[esys[[2,1]]]===List,1,0];
nevals=Length[esys[[2]]];
Join[esys[[1;;len,Complement[Range[nevals],remove]]],
If[Head[esys[[2,1]]]===List,{esys[[-1]]},{}]]
]



ClearAll[fix2PDmat];
fix2PDmat[mat_List?MatrixQ,\[Epsilon]_:10^-10]:=mat+DiagonalMatrix[ConstantArray[\[Epsilon],Length[mat]]]



ClearAll[PMapThread];
(*Parallel version of Mapthread that actually works*)
PMapThread[f_,exprs_List,lev_:1]:=WaitAll[MapThread[Composition[ParallelSubmit[Flatten[exprs],#]&,f],exprs,lev]]


(* ::Subsection:: *)
(*Main user functions*)


(* ::Subsubsection:: *)
(*EigenNDSolve*)


$dcon=DistributedContexts -> {"NMStools`Private`"};


Options[EigenNDSolve]={Npoly->80,returnInterp->True,noRowReduce->False,precision->MachinePrecision,
	returnMatrix->False,parallel->True};
EigenNDSolve::eqns="Supplied equations not in the correct form!";
EigenNDSolve::inhomog="Upgrade needed to deal with inhomogenous equations!";
EigenNDSolve::numeq="Number of variables does not match number of equations!";
EigenNDSolve::largeEl="Some derivative matrix elements are very large, `1`, and may be ill represented by machine precision numbers. Could increase Option, precision.";
EigenNDSolve::bounds="Unexpected boundary conditions supplied: `1`";
EigenNDSolve::noVBC="Not performing row reduction on boundary conditions due to unexpected boundary
	condition properties. Spurious eigenvalues will probably be produced";
EigenNDSolve::findzeros="Unable to find correct combination of non-zero elements to allow 
row-reduction. Proceeding calculation with full matrix";

EigenNDSolve[eqns_List,vars_List,range_List,evalvar_,OptionsPattern[]]:=Module[{z(*z is internal independent variable to use*),eqnBC,eqn,BC,dord,varQ,
eqnmat(*Terms multiplying each derivative of each variable*),BCnum,rowRedQ,inhomog,
T,dt,brules(*Abstract Chebyshev polynomials*),createOpMat,Npol=OptionValue[Npoly],Nprec=OptionValue[precision],sizelist,
OPnBC,BCrows,lv,chebvs,OPfull,
bcolind,browind,zeropos,choice,rowRinds,
BCfinal,BCSolfunc,OPfinal,partsfun,
eigsys,
pMap=If[OptionValue[parallel],ParallelMap[#1,#2,$dcon]&,Map[#1,#2]&]},

If[(*so it isn't necessary to find the matrix again in parameter scans*)
OptionValue[returnMatrix]&&MatrixQ@eqns[[1]],OPfinal=eqns[[1]];rowRedQ=False;lv=eqns[[2]],
(*Normal operation*)
(*Put on a {-1,1} region*)
varQ=(MemberQ[vars,#]&);(*Test is a variable is a var*)
eqnBC=transCo[eqns,varQ,range,{z,-1,1}];
(*Split into equations and boundary conditions*)
eqn=Select[eqnBC,MemberQ[#,f_?varQ[z]|Derivative[_][f_?varQ][z],Infinity]&];
BC=Complement[eqnBC,eqn];

If[Length[eqn]!= Length[vars],Message[EigenNDSolve::numeq];Abort[]];
(*Maximum derivative order list for each variable*)
dord=Max[Cases[eqns,Derivative[n_][#][_]:>n,Infinity]]&/@vars;
(*Convert eqns into single expression that equals zero*)
{eqn,BC}=If[And@@(Head[#]===Equal&/@#),Expand/@(Subtract@@@#),Message[EigenNDSolve::eqns];Abort[];]&/@{eqn,BC};
(*Check boundary conditions and calculate number of extra polynomials for each variable*)
{BCnum,rowRedQ}=checkBCform[BC,vars,dord];
If[OptionValue[noRowReduce],rowRedQ=False;];

(*Put equations into "matrix" form*)
eqnmat=Function[{expr},Function[{var},Coefficient[expr,Derivative[#][var][z]]&/@Range[0,Max@@dord]]/@vars]/@eqn;
inhomog=Join[eqn,BC]/.{f_?varQ[_]->0,Derivative[_][f_?varQ][_]->0};
If[Or@@(0=!=#&/@Chop[inhomog]),Message[EigenNDSolve::inhomog];Print[inhomog];Abort[]];

(*Define abstract Chebyshev properties to use in constructing matrix*)
brules=defineAbCheb[T,dt];

(*Function to make matrices*)
createOpMat[coef_,size_List]:=Block[{do=Length[coef],consFromT,nozmats,zfuncs,chebexpans,phivec,diffeq},
	(*Size is {N,N+BC} so actual matix size will be {N+1,N+BC+1}*)
	consFromT={Complement[Range[do],#],#}&@Select[Range[do],MemberQ[coef[[#]],z,{0,Infinity}]&];
	(*Calculate mat for those that don't contain z, this is seperate for speed*)
	nozmats=Total[coef[[#]]constructDmat[#-1,size+1]&/@consFromT[[1]]];
	If[nozmats===0,nozmats=SparseArray[{1,1}->0,size+1]];
	If[Max@Abs[Flatten[nozmats]]>10^17&&OptionValue[precision]===MachinePrecision,
			Message[EigenNDSolve::largeEl,Max@Abs[Flatten[nozmats]]>10^17]];

	(*Calculate for those that do contain z*)
	(*First, get Chebyshev expansions*)
	zfuncs=Function[Evaluate@{z},Evaluate[#]]&/@coef[[consFromT[[2]]]];
	phivec=(T[#]&/@Range[0,size[[2]]]);
	chebexpans=Map[Dot[#,phivec]&,Chop[CCapprox[#,size[[2]],evalvar,Nprec],10^-14]&/@zfuncs];
	(*Plug this into Chebyshev expansion*)
	diffeq=Total[pMap[ExpandAll,
		chebexpans (Nest[pMap[Expand[dt[#]]&,#]&,phivec,#]&/@(consFromT[[2]]-1))]];
	nozmats+Join[{diffeq/.T[_]->0},
		pMap[Coefficient[diffeq,T[#]]&, Range[size[[1]]]]]
];

(*Create full matrix without Boundary conditions*)
sizelist=ConstantArray[{Npol,#}&/@(Npol+BCnum),Length[BCnum]];
(*Print[Parallelize@Map[(Pause[0.2];#)&,Range[16]]//AbsoluteTiming];*)

(*This is just an implementation of MapThread[f,l,2] using only Map, as Parallelizing was somewhat broken*)
OPnBC=N[Chop[ArrayFlatten[Map[createOpMat@@##&,Transpose[{eqnmat,sizelist},{3,1,2}],{2}(*,
		DistributedContexts -> {"NMStools`Private`"}*)]],10^-14],Nprec];


(*Formulate boundary conditions and add into matrix*)
lv=Npol+BCnum;
chebvs=Map[Flatten@MapThread[#1(T[#]&/@Range[0,#2])&,{vars,lv}]/.#->1/.Thread[vars->0]&,vars];
BCrows=(BC/.Thread[vars->chebvs])/.{f_List[pm1_]:>(f/.brules[[-pm1]]),Derivative[n_][f_List][pm1_]:>(Nest[Map[Expand[dt[#]]&,#]&,f,n]/.brules[[-pm1]])};
OPfull=Join[OPnBC,N[BCrows,Nprec]];

If[rowRedQ,
(*First need to find {row,column} boundary points to be able to do a reduction*)
(*Indices of "tau" polynomials*)
bcolind=Sort@Flatten@MapThread[#1-Range[#2]&,{Rest@Accumulate[Prepend[lv+1,1]],BCnum}];
browind=Length[vars](Npol+1)+Range[Total[BCnum]];
(*Positions of zeros in each column*)
zeropos=Flatten@Position[OPfull[[browind,#]],a_/;a!=0]&/@bcolind;
(*Choose a combination that gives all 4 numbers without zero*)
choice=Flatten@Select[Tuples[zeropos],Sort[#]==Range[Length[zeropos]]&,1];


If[Length@choice==0,
(*In case it doesn't work, keep the program working*)
Message[EigenNDSolve::findzeros];OPfinal=OPfull,

rowRinds={browind[[choice]],bcolind}\[Transpose];
OPfull=Fold[RowR,OPfull,rowRinds];

BCfinal=OPfull[[browind]];
(*Function to solve for BCs given solution vector*)
BCSolfunc[solvec_]:=Block[{aBC,nbc=Range@Length[bcolind],fullsolvec},
fullsolvec=Fold[Insert[#1,#2[[2]],#2[[1]]]&,solvec,{bcolind,aBC/@nbc}\[Transpose]];
fullsolvec/.First@Solve[Chop[BCfinal.{fullsolvec}\[Transpose]]==0,aBC/@nbc]
];

partsfun=Complement[Range[Total[lv+1]],#]&;
OPfinal=OPfull[[partsfun[rowRinds[[All,1]]],partsfun[rowRinds[[All,2]]]]]
],

OPfinal=OPfull
];
];

If[(*Return matrix for parameter scans etc.*)
	OptionValue[returnMatrix]&&!MatrixQ[eqns[[1]]],{OPfinal,lv},
(*Normal operation*)
eigsys=Eigensystem[{-OPfinal/.evalvar->0,Coefficient[OPfinal,evalvar]}];
If[rowRedQ,
eigsys[[2]]=Map[BCSolfunc,eigsys[[2]]]];
eigsys[[2]]=Internal`PartitionRagged[#,lv+1]&/@eigsys[[2]];
(*Order according to size of imaginary part*)
eigsys=Join[{eigsys[[1]]},eigsys[[2]]\[Transpose]][[All,#]]&@
	Ordering[eigsys[[1]],All,
		If[Head[#1]===Complex&&Head[#2]===Complex,Im[#1]>Im[#2],Abs[#2]>Abs[#1]]&];

If[OptionValue[returnInterp],
eigsToInterp[eigsys,Rest@range],
eigsys]
]
]



(* ::Subsubsection:: *)
(*pseudoModes*)


(*pseudo modes returns functions to calculate the spatial structure and growth of the fastest growing pseudo-modes at a given time*)
(*esys is output of EigenNDSolve (with returnInterp\[Rule]False)*)
(*inPfunc is the inner product for energy, inP[{{u1,u2},{v1,v2},{w1,w2},...] for variables {u,v,w,...}*)
(*Output functions,
G[t] = maximal growth,
\[Xi][t] = perturbation that enables maximal growth, takes two time arguments (the second is the time at which Subscript[\[Kappa], 0] is evaluated)*)
pseudoModes::nointerp="Function requires returnInterp\[Rule]False in EigenNDSolve!";
Options[pseudoModes]={removeLarEvalsCutoff->10^3,HermitianFix\[Epsilon]->0,cholPrecision->MachinePrecision,returnICs->False,
	returnRedA->False,parallel->True};
pseudoModes[esys_List,Kkeep_Integer,inPfunc_,OptionsPattern[]]:=Module[{topKmodes,normmat,FUmat,Mmat,Mmat1,
Fmat,iFmat,e\[CapitalLambda]t,\[Kappa]0fun,
pMap=If[OptionValue[parallel],ParallelMap[#1,#2,$dcon]&,Map[#1,#2]&]},
If[Head[esys[[2,1]]]=!=List,Message[pseudoModes::nointerp]];
(*Remove unphysical eigenmodes and choose top K modes*)
topKmodes=removeLargeEvalues[esys,OptionValue[removeLarEvalsCutoff]];
topKmodes=topKmodes[[All,
	Ordering[Im[topKmodes[[1]]],Kkeep,#1>#2&]]];

(*Calculate matrix of inner products*)
If[TrueQ@MatrixQ[inPfunc],
(*If a matrix is supplied, assume this is inner product matrix for all variables*)
	normmat=CholeskyDecomposition[inPfunc];
	FUmat=ArrayFlatten[{pMap[normmat.{#}\[Transpose]&,
		(Join@@topKmodes[[2;;,#]])&/@Range[Length[topKmodes[[2]]]]]}];
	FUmat=SetPrecision[FUmat,OptionValue[cholPrecision]];
	Mmat=FUmat\[ConjugateTranspose].FUmat,

	topKmodes=SetPrecision[topKmodes,OptionValue[cholPrecision]];
	Mmat=Chop[Table[
		inPfunc[{#[[i]],#[[j]]}&/@Rest[topKmodes]],
		{i,1,Kkeep},{j,1,Kkeep}]]
];
Fmat=CholeskyDecomposition[fix2PDmat[Mmat,OptionValue[HermitianFix\[Epsilon]]]];
iFmat=Inverse[Fmat];
(*Set back to MachinePrecision to make G and \[Xi] faster to evaluate*)
{Fmat,iFmat}=SetPrecision[{Fmat,iFmat},MachinePrecision];
If[OptionValue[returnRedA]===True,Fmat.DiagonalMatrix[-I topKmodes[[1]]].iFmat,
	e\[CapitalLambda]t=Function[{t},DiagonalMatrix[Exp[-I topKmodes[[1]] t]]];

	(*{G[t],\[Xi][t]*)
	\[Kappa]0fun=Function[{t},Flatten[iFmat.(SingularValueDecomposition[Fmat.e\[CapitalLambda]t[t].iFmat,1][[3]])]];
	{Function[{t},SingularValueList[Fmat.e\[CapitalLambda]t[t].iFmat,1][[1]]^2],
		Function[{t,t0},(Flatten[(e\[CapitalLambda]t[t].{\[Kappa]0fun[t0]}\[Transpose])].topKmodes[[#]])&/@Range[2,Length[topKmodes]]],
		If[OptionValue[returnICs],\[Kappa]0fun]/.Null->Sequence[]}
]

]


(* ::Subsubsection::Closed:: *)
(*Plotting*)


(*Plotting e-values*)
plotnum=1;
Options[plotEValues]={estateList->False};
plotEValues::estateform="E-states must be given in the form of an interpolating function";
plotEValues[evallist_,range_,OptionsPattern[]]:=DynamicModule[{imre,funrange},
imre={Re[#],Im[#]}&;

If[OptionValue[estateList]=!=False,
If[Head[OptionValue[estateList][[1,1]]]=!=InterpolatingFunction,Message[plotEValues::estateform],
funrange=OptionValue[estateList][[1,1,1]]]];

Show[
ListPlot[MapIndexed[Button[Tooltip[#1,ToString[#2[[1]]]<>": "<>ToString[#1]],plotnum=#2[[1]]]&,imre/@evallist],PlotStyle->Directive[PointSize[Medium],Black],PlotRange->range],
Epilog->Inset[  
If[OptionValue[estateList]=!=False,
Dynamic[
intFuncsPlot[OptionValue[estateList][[All,plotnum]]]
],
Graphics[{White,Point[{0,0}]}]],
ImageScaled[{0,0}],ImageScaled[{0,0}],0.5(range[[1,2]]-range[[1,1]])]
]


]


(* ::Section::Closed:: *)
(*Old functions*)


(*Returns eigensystem of a generalized polynomial eigenvalue problem*)
Options[polyEigen]={computeEstates->True};
polyEigen::unmatched="Expected order does not match with that found using matrix! Orders `1`";
polyEigen::matsize="F or G have come out to the wrong size for some reason";
polyEigen[eqnmat_,{ename_,ordeval_Integer},OptionsPattern[]]:=Module[{matlist,zerocheck,indsiz,Fmat,Gmat,esys},
(*Break up matrix into \[Omega] polynomial*)
matlist=Coefficient[eqnmat,ename,#]&/@Range[0,ordeval]//N;
(*Check if any of the coefficients are zero*)
zerocheck=Function[{arr},And@@(0==#&/@ArrayRules[arr][[All,2]])]/@matlist;
If[Or@@zerocheck,Message[polyEigen::unmatched,zerocheck]];

(*Form larger matrices*)
indsiz=Length[matlist[[1]]];
If[ordeval==1,
Fmat=matlist[[1]];
Gmat=-matlist[[2]];,

Fmat=SparseArray`SparseBlockMatrix[Join[{Band[{1,2},{ordeval-1,ordeval}]->IdentityMatrix[indsiz]},({ordeval,#}->-matlist[[#]]&/@Range[ordeval])],{ordeval,ordeval}];
Gmat=SparseArray[Band[{1,1}]->Join[ConstantArray[IdentityMatrix[indsiz],ordeval-1],{matlist[[ordeval+1]]}]];
];

If[(Dimensions/@{Fmat,Gmat})!={indsiz ordeval,indsiz ordeval},Message[polyEigen::matsize]];

If[!OptionValue[computeEstates],
esys=Eigenvalues[{Fmat,Gmat}];Sort[esys],
(*Only want the first indsiz terms of each eigenvector*)
esys=Eigensystem[{Fmat,Gmat}];
(*Pick off first bit of each estate and re-order*)
Function[{ord},{esys[[1,ord]],(Take[#,-indsiz]&/@esys[[2]])[[ord]]}]@Ordering[Abs@esys[[1]]]
]
]



Options[FDEigenNDSolve]={gridN->201,diffOrder->4,computeEstates->True,returnInterp->True,returnFDMatrix->False};
FDEigenNDSolve::eqns="Supplied equations not in the correct form!";
FDEigenNDSolve::inhomog="Upgrade needed to deal with inhomogenous equations!";
FDEigenNDSolve::numeq="Number of variables does not match number of equations!";
FDEigenNDSolve[eqns_List,vars_List,range_List,evalvar_List,OptionsPattern[]]:=
Module[{eqn0,eqnmat,inhomog,dord,
x,np,r1,r2,\[Delta]r,grid,ingrid,consvec,spD,coeffmats,
fullop,esystem},

x=range[[1]];
If[Length[eqns]!= Length[vars],Message[FDEigenNDSolve::numeq];Abort[]];
(*Maximum derivative order*)
dord=Max@Cases[eqns,Derivative[n_][_][_]:>n,Infinity];

(*Convert eqns into single expression that equals zero*)
eqn0=If[And@@(Head[#]===Equal&/@eqns),Expand/@(Subtract@@@eqns),Message[FDEigenNDSolve::eqns];Abort[];];

(*Put equations into "matrix" form*)
eqnmat=Function[{expr},Function[{var},Coefficient[expr,Derivative[#][var][x]]&/@Range[0,dord]]/@vars]/@eqn0;
inhomog=Expand[
eqn0-Plus@@@
(Function[{dlist},
MapThread[Function[{obj,var},Plus@@(obj[[#+1]]Derivative[#][var][x]&/@Range[0,dord])],{dlist,vars}]]/@
eqnmat)
];
If[Chop[Plus@@inhomog]=!=0,Message[FDEigenNDSolve::inhomog];Print[inhomog]];

(*Setup grid*)
r1=range[[2]];r2=range[[3]];
np=OptionValue[gridN];\[Delta]r=(r2-r1)/(np-1);
grid=  N@Range[r1,r2,\[Delta]r];ingrid=remEnds@grid;

(*Function to create a finite difference matrix for given derivative order*)
dmatF=CoefficientArrays[remEnds@(NDSolve`FiniteDifferenceDerivative[Derivative[#],grid, Map[afun, grid],"DifferenceOrder"->OptionValue[diffOrder]]/.{afun[N@r1]->0,afun[N@r2]->0}),
afun/@ingrid][[2]]&;
(*Create base matrices to construct full differential operator*)
consvec=ConstantArray[#,{np-2}]&;
spD=SparseArray[Band[{1,1}]->#]&;

(*Construct full operator*)
coeffmats=Map[spD,eqnmat/.{evalvar[[1]]->consvec[evalvar[[1]]],x->ingrid,num_?NumericQ:>consvec[num]},{3}];

fullop=
ArrayFlatten[
Map[
Function[{coef},Plus@@MapThread[#1.#2&,{coef,Join[{spD@consvec@1},dmatF/@Range[dord]]}]] ,
coeffmats,{2}]
];

(*with returnFDMatrix\[Rule]True, no eigenvalues are solved for and it just returns the finite difference matrix*)
If[OptionValue[returnFDMatrix]===False,

(*Find eigenstates using polyEigen*)
esystem=polyEigen[fullop,evalvar,computeEstates->OptionValue[computeEstates]];
(*Split up eigenstates according to variables*)
If[OptionValue[computeEstates]===True,

esystem=Join[{esystem[[1]]},Map[addEnds,Transpose[Partition[#,np-2]&/@esystem[[2]]],{2}],{grid}];

If[TrueQ@OptionValue[returnInterp],
Join[{esystem[[1]]},Map[ListInterpolation[#,{esystem[[-1]]}]&,remEnds@esystem,{2}]],
esystem],


esystem],

fullop]



]




(*DOESN'T WORK ON ANY OF THE SYSTEMS I'M LOOKING AT SINCE THEY ARE DIFFERENTIAL ALGEBRAIC PDEs!!!*)
Options[linearNPDESolve]={gridN->201,diffOrder->4};
(*linearNPDESolve solves the linear PDE to compare with eigenvalue solutions*)
(*eqns must be given in form \!\(
\*SubscriptBox[\(\[PartialD]\), \(t\)]\[Xi]\)=\[ScriptCapitalL]\[Xi], where \[ScriptCapitalL]\[Xi] is supplied as the first argument*)
(*vars and xrange as for FDEigenNDSolve*)
(*trange is a list containing {tmax,\[Delta]t}*)
(*IC is a list of the initial conditions*)
linearNPDESolve::initL="Size of initial conditions does not match gridN";
linearNPDESolve[eqns_List,vars_List,xrange_List,trange_List,IC_List,OptionsPattern[]]:=Module[
{\[Delta]t=trange[[2]],tmax=trange[[1]],eqnsin=Thread[eqns==0](*So as to easily use FDEigenNDSolve*),localevar(*Just a dummy*),initC,
FDmat,RKts,fullsol},

If[Length[IC[[1]]]=!=OptionValue[gridN],Message[linearNPDESolve::initL];Abort,initC=Flatten[remEnds/@IC]];

(*Find the finite difference matrix*)
FDmat=FDEigenNDSolve[eqnsin,vars,xrange,{localevar,1},returnFDMatrix->True,gridN->OptionValue[gridN],diffOrder->OptionValue[diffOrder]];


RKts[\[Xi]n_]:=Block[{k1,k2,k3,k4},
k1=FDmat.\[Xi]n;
k2=FDmat.(\[Xi]n+\[Delta]t/2 k1);
k3=FDmat.(\[Xi]n+\[Delta]t/2 k2);
k4=FDmat.(\[Xi]n+\[Delta]t k3);
\[Xi]n+\[Delta]t/6 (k1+2k2+2k3+k4)
];


fullsol=NestList[RKts,{initC}\[Transpose],Floor[tmax/\[Delta]t]]//myTiming;
Print[Dimensions[Map[addEnds,
(Partition[#,OptionValue[gridN]-2]&/@fullsol)\[Transpose],{2}]]];

ListInterpolation[#,{{0,tmax},Rest@xrange}]&/@
(Map[addEnds,
(Partition[#,OptionValue[gridN]-2]&/@fullsol)\[Transpose],{2}])//myTiming

]



End[]
EndPackage[]
