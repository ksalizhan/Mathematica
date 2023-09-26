(* ::Package:: *)

BeginPackage["defs448`"]

conj::usage="Replaces I by -I in an expression, also \!\(\*SuperscriptBox[\(a\), \(*\)]\)"
(*\[CircleTimes]::"Kronecker or Outer Product"*)
MF::usage="Better Display of numeric Matrices"
ThadPlot::usage="Plot[a,b,Style options,{xlabel,ylabel}] or Show[oldplot,{xlabel,ylabel}]"
(* a^\[Dagger] is the Hermitian conjugate (conjugate transpose)*)
norm::usage="normalize a vector"
\[DoubleStruckOne]::usage="Identity Matrix \!\(\*SubscriptBox[\"II\", \"g\"]\)"
angmom::usage="generate {\!\(\*SubscriptBox[\(J\), \(x\)]\),\!\(\*SubscriptBox[\(J\), \(y\)]\),\!\(\*SubscriptBox[\(J\), \(z\)]\)}"
P::usage="Hydrogen wave functions, atomic units \!\(\*SubscriptBox[\(P\), \(n, l\)]\)[r]"
\[Epsilon]::usage="totally antisymmetric tensor in 3d"
\[Delta]::usage="Kronecker delta"
Begin["`Private`"]
<<Notation`
conj[x_]:=x/.Complex[a_,b_]->Complex[a,-b]
MF[q_]:=Quiet[NumberForm[MatrixForm[Chop[q]],3]]
ThadPlot[a_,b_,c_,d_]:=Plot[a,b,{c,AxesStyle->Directive[Black,18],GridLines->Automatic,GridLinesStyle->LightGray,Frame->True,Axes->False,
	FrameLabel->{Text[Style[d[[1]],FontSize->14]],Text[Style[d[[2]],FontSize->14]]},FrameStyle->{{Directive[GrayLevel[0],14],GrayLevel[1]},{Directive[GrayLevel[0],14],GrayLevel[1]}}}]
ThadPlot[a_,b_,c_]:=ThadPlot[a,b,PlotStyle->Thick,c];
ThadPlot[a_,b_]:=Show[a,{AxesStyle->Directive[Black,18],GridLines->Automatic,GridLinesStyle->Directive[Thin,LightGray],Axes->False,Frame->True,FrameLabel->b,FrameStyle->{{Directive[GrayLevel[0],14],GrayLevel[1]},{Directive[GrayLevel[0],14],GrayLevel[1]}}}]
a_\[CircleTimes]b_:=If[Length[Dimensions[a]]==Length[Dimensions[b]]==1,Flatten[KroneckerProduct[{a},{b}]], KroneckerProduct[a,b]]
a_\[CircleTimes]b_\[CircleTimes]c_:=(a\[CircleTimes]b)\[CircleTimes]c
SuperStar[a_]:=a/.Complex[x_,y_]->Complex[x,-y]
SuperDagger[b_]:=Transpose[SuperStar[b]]
norm[v_]:=v/Sqrt[v.SuperStar[v]]
angmom[j_]:=With[{jp=DiagonalMatrix[Sqrt[j(j+1)-ms(ms+1)]/.ms->Range[j-1,-j,-1],1]},
{(jp+SuperDagger[jp])/2,(jp-SuperDagger[jp])/(2 I),DiagonalMatrix[Range[j,-j,-1]]}
]
Ket[a_]:={norm[Flatten[a]]}\[Transpose]//Simplify
Bra[a_]:=SuperDagger[Ket[a]]//Simplify
BraKet[a_,b_]:=Bra[Flatten[a]].Ket[Flatten[b]]//Tr
Subscript[\[DoubleStruckOne], n_]:=IdentityMatrix[n];
Subscript[\[Epsilon], i_,j_,k_]:=Signature[{i,j,k}];
Subscript[\[Delta], i_,j_]:=KroneckerDelta[i,j];
w[\[Kappa]_,\[Mu]_,z_]:=Exp[-z/2]z^(\[Mu]+1/2) HypergeometricU[1/2+\[Mu]-\[Kappa],1+2\[Mu],z];
Subscript[P, \[Nu]_,l_][r_]:=1/Sqrt[\[Nu]^2 Gamma[\[Nu]+l+1]Gamma[\[Nu]-l]] w[\[Nu],l+1/2,2 r/\[Nu]];
End[]
EndPackage[]

















