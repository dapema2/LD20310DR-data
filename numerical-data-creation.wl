(* ::Package:: *)

(* ::Section:: *)
(*Parameters*)


foldername = ""; (* The name of the subfolder where the data will be saved *)
filename = ""; (* The name of the files where the data will be stored (DON'T WRITE FILE EXTENSIONS) *)

prec = 155; (* precision (high required for large a) *)
cL = 0.5;
cS =1;
cR =2;
v = -cS; 
m = 1;
\[HBar] = 1;
a = {0, 1, 10, 20, 50, 100}; (* values of central width *)

wmin = 10^-4;
wmax = wmaximum;
wnum = 10000;

(* for constant separation of log(w) *)
wloglist = Table[N[logw, prec], {logw, Log[10, wmin], Log[10, wmax], (Log[10, wmax] - Log[10,wmin])/(wnum-1)}];
wlist  = Table[N[10^(wloglist[[i]]), prec], {i, Length[wloglist]}];

(* Don't touch anything below this line *)
\[Xi]L = \[HBar]/(m cL); \[Xi]S = \[HBar]/(m cS); \[Xi]R = \[HBar]/(m cR);
kmax = -1/\[Xi]L \[Sqrt](-2+v^2/(2 cL^2) + Abs[v]/(2 cL) Sqrt[8+v^2/cL^2]);
wmaximum = v kmax - cL \[Sqrt](kmax^2 + (\[Xi]L^2 kmax^4)/4)


(* ::Section:: *)
(*Wavevectors and Coefficients*)


(* Create the equations *)
subeq = Table[(w - v k)^2 == cR^2 (k^2+(k^4 \[Xi]R^2)/4), {w, wlist}];
supereq = Table[(w - v k)^2 == cL^2 (k^2+(k^4 \[Xi]L^2)/4), {w, wlist}];
soniceq = Table[(w + cS k)^2 == cS^2 (k^2+(k^4 \[Xi]S^2)/4), {w, wlist}];

(* Solve the equations *)
ksub = Table[k /. NSolve[subeq[[i]], k, WorkingPrecision->prec], {i, Length[subeq]}];
ksuper =  Table[k /. NSolve[supereq[[i]], k, WorkingPrecision->prec], {i, Length[supereq]}];
ksonic =  Table[k /. NSolve[soniceq[[i]], k, WorkingPrecision->prec], {i, Length[soniceq]}];

(* Rearrange the solutions in the order (v,u,+,-) *)
ksub = Table[SortBy[ksub[[j]],{Boole[Im[#]!=0],(*Real numbers first (Im==0)*)-Im[#](*Complex numbers sorted by decreasing imaginary part*)}&], {j, Length[wlist]}];
ksuper = Table[SortBy[ksuper[[j]],{Boole[Im[#]!=0],(*Real numbers first (Im==0)*)-Im[#](*Complex numbers sorted by decreasing imaginary part*)}&], {j, Length[wlist]}];
ksuper = ksuper[[All, {3,2,4,1}]];
ksonic = Table[SortBy[ksonic[[j]],{Boole[Im[#]!=0],(*Real numbers first (Im==0)*)-Im[#](*Complex numbers sorted by decreasing imaginary part*)}&], {j, Length[wlist]}];

(* Derivatives of the wavevectors *) 
kderivative = (2 (w - v k))/(c^2 (2 k + k^3 \[Xi]^2) + 2v(w - v k)); (* Found by derivating the equations *)
dksub = Table[kderivative /. {k -> ksub[[i]], c -> cR, \[Xi] -> \[Xi]R, w -> wlist[[i]]}, {i, Length[wlist]}];
dksuper= Table[kderivative /. {k -> ksuper[[i]], c -> cL, \[Xi] -> \[Xi]L, w -> wlist[[i]]}, {i, Length[wlist]}];
dksonic = Table[kderivative /. {k -> ksonic[[i]], c -> cS, \[Xi] -> \[Xi]S, w -> wlist[[i]]}, {i, Length[wlist]}];

(* Normalisation Coefficients *)
n = 1/(4 \[Pi] \[HBar]);
Dcoef=((w - v k + (c \[Xi] k^2)/2)Sqrt[Abs[dk]])/(\[Sqrt](4 \[Pi] \[HBar] n \[Xi] c Abs[w (Re[k]^2-Im[k]^2)-v Re[k] (Re[k]^2+Im[k]^2)]));
Ecoef=-(((w - v k - (c \[Xi] k^2)/2)Sqrt[Abs[dk]])/(\[Sqrt](4 \[Pi] \[HBar] n \[Xi] c Abs[w (Re[k]^2-Im[k]^2)-v Re[k] (Re[k]^2+Im[k]^2)])));
Dsub = Table[Dcoef /. {k -> ksub[[i]], dk -> dksub[[i]], c -> cR, \[Xi] -> \[Xi]R, w -> wlist[[i]]}, {i, Length[wlist]}];
Esub = Table[Ecoef /. {k -> ksub[[i]], dk -> dksub[[i]], c -> cR, \[Xi] -> \[Xi]R, w -> wlist[[i]]}, {i, Length[wlist]}];
Dsuper = Table[Dcoef /. {k -> ksuper[[i]], dk -> dksuper[[i]], c -> cL, \[Xi] -> \[Xi]L, w -> wlist[[i]]}, {i, Length[wlist]}];
Esuper = Table[Ecoef /. {k -> ksuper[[i]], dk -> dksuper[[i]], c -> cL, \[Xi] -> \[Xi]L, w -> wlist[[i]]}, {i, Length[wlist]}];
Dsonic = Table[Dcoef /. {k -> ksonic[[i]], dk -> dksonic[[i]], c -> cS, \[Xi] -> \[Xi]S, w -> wlist[[i]]}, {i, Length[wlist]}];
Esonic = Table[Ecoef /. {k -> ksonic[[i]], dk -> dksonic[[i]], c -> cS, \[Xi] -> \[Xi]S, w -> wlist[[i]]}, {i, Length[wlist]}];

(* Clear garbage variables *)
Clear[{subeq, supereq, soniceq, ksubaux, ksuperaux, ksonicaux, kderivative, dksub, dksuper, dksonic, n, Dcoef, Ecoef}];


(* ::Subsection:: *)
(*Verification of order (manual)*)


(* In this cell we check if the equations are correctly ordered using the low-w expansions, just in case. The canonical order must be {v,u,+,-}, where (+ -> 3) and (- -> 4) for the supersonic region *)
(* VERY IMPORTANT! The script will fail if the order is not correct! *)
x = 1;
Print["w = ", N[wlist[[x]]]];
Print[];
(* subsonic *)
Print["SUBSONIC"];
Print["Current Order: ", N[ksub[[x]]]]
Print["Correct Order: ", N[Normal[{SeriesData[w, 0, {(-c + v)^(-1), 0, Rational[1, 8] c (c - v)^(-4) \[Xi]^2}, 1, 4, 1],SeriesData[w, 0, {(c + v)^(-1), 0, Rational[-1, 8] c (c + v)^(-4) \[Xi]^2}, 1, 4, 1],SeriesData[w, 0, {Complex[0, 2] c^(-1) (((c - v) (c + v))^Rational[1, 2]/\[Xi]), v/(c^2 - v^2), Complex[0, Rational[1, 4]] c ((c - v) (c + v))^Rational[-5, 2] (c^2 + 2 v^2) \[Xi], Rational[-1, 2] c^2 v (c^2 - v^2)^(-4) (c^2 + v^2) \[Xi]^2}, 0, 4, 1],SeriesData[w, 0, {Complex[0, -2] c^(-1) (((c - v) (c + v))^Rational[1, 2]/\[Xi]), v/(c^2 - v^2), Complex[0, Rational[-1, 4]] c ((c - v) (c + v))^Rational[-5, 2] (c^2 + 2 v^2) \[Xi], Rational[-1, 2] c^2 v (c^2 - v^2)^(-4) (c^2 + v^2) \[Xi]^2}, 0, 4, 1]}]/. {c -> cR, \[Xi] -> \[Xi]R, w -> wlist[[x]]}]];
Print[];
(* supersonic *)
Print["SUPERSONIC"];
Print["Current Order: ", N[ksuper[[x]]]];
Print["Correct Order: ", N[Normal[{SeriesData[w, 0, {(-c + v)^(-1), 0, Rational[1, 8] c (c - v)^(-4) \[Xi]^2}, 1, 4, 1],SeriesData[w, 0, {(c + v)^(-1), 0, Rational[-1, 8] c (c + v)^(-4) \[Xi]^2}, 1, 4, 1],SeriesData[w, 0, {2 c^(-1) ((-c^2 + v^2)^Rational[1, 2]/\[Xi]), v/(c^2 - v^2), Rational[-1, 4] c (-c^2 + v^2)^Rational[-5, 2] (c^2 + 2 v^2) \[Xi], Rational[-1, 2] c^2 v (c^2 - v^2)^(-4) (c^2 + v^2) \[Xi]^2}, 0, 4, 1],SeriesData[w, 0, {(-2) c^(-1) ((-c^2 + v^2)^Rational[1, 2]/\[Xi]), v/(c^2 - v^2), Rational[1, 4] c (-c^2 + v^2)^Rational[-5, 2] (c^2 + 2 v^2) \[Xi], Rational[-1, 2] c^2 v (c^2 - v^2)^(-4) (c^2 + v^2) \[Xi]^2}, 0, 4, 1]} ]/. {c -> cL, \[Xi] -> \[Xi]L, w -> wlist[[x]]}]];
Print[];
(* sonic *)
Print["SONIC"];
Print["Current Order: ", N[ksonic[[x]]]];
Print["Correct Order: ", N[Normal[{SeriesData[w, 0, {Rational[-1, 2]/c, 0, 0, 0, 0, 0, Rational[1, 384] c^(-3) \[Xi]^2}, 3, 10, 3],SeriesData[w, 0, {2 (c \[Xi]^2)^Rational[-1, 3], 0, Rational[1, 6]/c, 0, Rational[-1, 36] c^Rational[-5, 3] \[Xi]^Rational[2, 3], 0, Rational[5, 648] c^Rational[-7, 3] \[Xi]^Rational[4, 3], 0, Rational[-1, 384] c^(-3) \[Xi]^2}, 1, 10, 3],SeriesData[w, 0, {(-1 + Complex[0, 1] 3^Rational[1, 2]) (c \[Xi]^2)^Rational[-1, 3], 0, Rational[1, 6]/c, 0, Rational[1, 72] (1 + Complex[0, 1] 3^Rational[1, 2]) c^Rational[-5, 3] \[Xi]^Rational[2, 3], 0, Rational[5, 1296] (-1 + Complex[0, 1] 3^Rational[1, 2]) c^Rational[-7, 3] \[Xi]^Rational[4, 3], 0, Rational[-1, 384] c^(-3) \[Xi]^2}, 1, 10, 3],SeriesData[w, 0, {(-1 + Complex[0, -1] 3^Rational[1, 2]) (c \[Xi]^2)^Rational[-1, 3], 0, Rational[1, 6]/c, 0, Rational[1, 72] (1 + Complex[0, -1] 3^Rational[1, 2]) c^Rational[-5, 3] \[Xi]^Rational[2, 3], 0, Rational[5, 1296] (-1 + Complex[0, -1] 3^Rational[1, 2]) c^Rational[-7, 3] \[Xi]^Rational[4, 3], 0, Rational[1, 384] c^(-3) \[Xi]^2}, 1, 10, 3]}]/. {c -> cS, \[Xi] -> \[Xi]S, w -> wlist[[x]]}]];


(* ::Section:: *)
(*Matching Equations*)


(* Matching equations with x = 0 in LS and x = a in SR *) 
matchingEquations = Table[{AL1 Dsuper[[i,1]]+AL2 Dsuper[[i,2]]+AL3 Dsuper[[i,3]]+AL4 Dsuper[[i,4]]==AS1 Dsonic[[i,1]]+AS2 Dsonic[[i,2]]+AS3 Dsonic[[i,3]]+AS4 Dsonic[[i,4]],I AL1 Dsuper[[i,1]] ksuper[[i,1]]+I AL2 Dsuper[[i,2]] ksuper[[i,2]]+I AL3 Dsuper[[i,3]] ksuper[[i,3]]+I AL4 Dsuper[[i,4]] ksuper[[i,4]]==I AS1 Dsonic[[i,1]] ksonic[[i,1]]+I AS2 Dsonic[[i,2]] ksonic[[i,2]]+I AS3 Dsonic[[i,3]] ksonic[[i,3]]+I AS4 Dsonic[[i,4]] ksonic[[i,4]],AL1 Esuper[[i,1]]+AL2 Esuper[[i,2]]+AL3 Esuper[[i,3]]+AL4 Esuper[[i,4]]==AS1 Esonic[[i,1]]+AS2 Esonic[[i,2]]+AS3 Esonic[[i,3]]+AS4 Esonic[[i,4]],I AL1 Esuper[[i,1]] ksuper[[i,1]]+I AL2 Esuper[[i,2]] ksuper[[i,2]]+I AL3 Esuper[[i,3]] ksuper[[i,3]]+I AL4 Esuper[[i,4]] ksuper[[i,4]]==I AS1 Esonic[[i,1]] ksonic[[i,1]]+I AS2 Esonic[[i,2]] ksonic[[i,2]]+I AS3 Esonic[[i,3]] ksonic[[i,3]]+I AS4 Esonic[[i,4]] ksonic[[i,4]],AS1 Dsonic[[i,1]] E^(I a[[j]] ksonic[[i,1]])+AS2 Dsonic[[i,2]] E^(I a[[j]] ksonic[[i,2]])+AS3 Dsonic[[i,3]] E^(I a[[j]] ksonic[[i,3]])+AS4 Dsonic[[i,4]] E^(I a[[j]] ksonic[[i,4]])==AR1 Dsub[[i,1]] E^(I a[[j]] ksub[[i,1]])+AR2 Dsub[[i,2]] E^(I a[[j]] ksub[[i,2]])+AR3 Dsub[[i,3]] E^(I a[[j]] ksub[[i,3]])+AR4 Dsub[[i,4]] E^(I a[[j]] ksub[[i,4]]),I AS1 Dsonic[[i,1]] E^(I a[[j]] ksonic[[i,1]]) ksonic[[i,1]]+I AS2 Dsonic[[i,2]] E^(I a[[j]] ksonic[[i,2]]) ksonic[[i,2]]+I AS3 Dsonic[[i,3]] E^(I a[[j]] ksonic[[i,3]]) ksonic[[i,3]]+I AS4 Dsonic[[i,4]] E^(I a[[j]] ksonic[[i,4]]) ksonic[[i,4]]==I AR1 Dsub[[i,1]] E^(I a[[j]] ksub[[i,1]]) ksub[[i,1]]+I AR2 Dsub[[i,2]] E^(I a[[j]] ksub[[i,2]]) ksub[[i,2]]+I AR3 Dsub[[i,3]] E^(I a[[j]] ksub[[i,3]]) ksub[[i,3]]+I AR4 Dsub[[i,4]] E^(I a[[j]] ksub[[i,4]]) ksub[[i,4]],AS1 E^(I a[[j]] ksonic[[i,1]]) Esonic[[i,1]]+AS2 E^(I a[[j]] ksonic[[i,2]]) Esonic[[i,2]]+AS3 E^(I a[[j]] ksonic[[i,3]]) Esonic[[i,3]]+AS4 E^(I a[[j]] ksonic[[i,4]]) Esonic[[i,4]]==AR1 E^(I a[[j]] ksub[[i,1]]) Esub[[i,1]]+AR2 E^(I a[[j]] ksub[[i,2]]) Esub[[i,2]]+AR3 E^(I a[[j]] ksub[[i,3]]) Esub[[i,3]]+AR4 E^(I a[[j]] ksub[[i,4]]) Esub[[i,4]],I AS1 E^(I a[[j]] ksonic[[i,1]]) Esonic[[i,1]] ksonic[[i,1]]+I AS2 E^(I a[[j]] ksonic[[i,2]]) Esonic[[i,2]] ksonic[[i,2]]+I AS3 E^(I a[[j]] ksonic[[i,3]]) Esonic[[i,3]] ksonic[[i,3]]+I AS4 E^(I a[[j]] ksonic[[i,4]]) Esonic[[i,4]] ksonic[[i,4]]==I AR1 E^(I a[[j]] ksub[[i,1]]) Esub[[i,1]] ksub[[i,1]]+I AR2 E^(I a[[j]] ksub[[i,2]]) Esub[[i,2]] ksub[[i,2]]+I AR3 E^(I a[[j]] ksub[[i,3]]) Esub[[i,3]] ksub[[i,3]]+I AR4 E^(I a[[j]] ksub[[i,4]]) Esub[[i,4]] ksub[[i,4]]}, {i, Length[wlist]}, {j, Length[a]}];


outputfiles = FileNameJoin[{NotebookDirectory[], foldername, filename}]

If[FileExistsQ[filename], DeleteFile[filename]];
Do[
scatteringMatrix4L = Table[ N[{AR2, AL1, AL2, AR3, AS1, AS2, AS3, AS4} /. Flatten[NSolve[matchingEquations[[\[Omega],\[Alpha]]] /. {AR4 -> 0, AR1 -> 0, AL3 -> 0, AL4 -> 1}, {AR2, AR3, AS1, AS2, AS3, AS4, AL1, AL2}, WorkingPrecision -> prec]],prec],{\[Omega], Length[wlist]}]; (* You can change the conditions {AR4 -> 0, AR1 -> 0, AL3 -> 0, AL4 -> 1} to change the in mode *)
output = FileNameJoin[{NotebookDirectory[], foldername, filename <> "a" <> ToString[\[Alpha]] <> ".txt"}];
If[FileExistsQ[output], DeleteFile[output]];
Export[outputTemp, scatteringMatrix4L, "Table"];
,{\[Alpha], Length[a]}
]
